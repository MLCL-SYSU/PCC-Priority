/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*  Special thanks to author ZhangSoonyang
 *  Our code is developed from the code repository https://github.com/SoonyangZhang/DrainQueueCongestion
 *  
 *  Authors:Zhuoyue Chen
 * 
 *  BSD 3-Clause License
 *  Copyright (c) 2021, SoonyangZhang
 *  All rights reserved.
 */
#include <algorithm>
#include "pcc_sender.h"
#include "unacked_packet_map.h"
#include "commandlineflags.h"
#include "rtt_stats.h"
#include "proto_time.h"
#include "random.h"
#include "quic_logging.h"

#define FALLTHROUGH_INTENDED \
  do {                            \
  } while (0)
namespace qcc {

DEFINE_bool(enable_rtt_deviation_based_early_termination, false,
            "Whether trigger early termination by comparing inflated rtt with "
            "rtt deviation");
DEFINE_bool(trigger_early_termination_based_on_interval_queue_front, false,
            "Whether trigger early termination by comparing most recent "
            "smoothed rtt and the rtt upon start of the interval at the front "
            "of the interval queue. The current interval in the queue is used "
            "by default");
DEFINE_bool(enable_early_termination_based_on_latest_rtt_trend, false,
            "Whether trigger early termination by comparing latest RTT and "
            "smoothed RTT");
DEFINE_double(max_rtt_fluctuation_tolerance_ratio_in_starting, 0.3,
              "Ignore RTT fluctuation within 30 percent in STARTING mode");
DEFINE_double(max_rtt_fluctuation_tolerance_ratio_in_decision_made, 0.05,
              "Ignore RTT fluctuation within 5 percent in DECISION_MADE mode");
DEFINE_double(rtt_fluctuation_tolerance_gain_in_starting, 2.5,
              "Ignore rtt fluctuation within 2.5 multiple of rtt deviation in "
              "STARTING mode.");
DEFINE_double(rtt_fluctuation_tolerance_gain_in_decision_made, 1.5,
              "Ignore rtt fluctuation within 1.5 multiple of rtt deviation in "
              "DECISION_MADE modes.");
DEFINE_bool(can_send_respect_congestion_window, false,
            "Use calculated congestion window to determine whether CanSend "
            "should return true.");
DEFINE_double(bytes_in_flight_gain, 2.5,
              "Enable a specific multiple of approxmate cwnd bytes in flight.");
DEFINE_bool(exit_starting_based_on_sampled_bandwidth, false,
            "When exiting STARTING, fall back to the minimum of the max "
            "bandwidth by bandwidth sampler and half of current sending rate");
DEFINE_bool(restore_central_rate_upon_app_limited, false,
            "Whether restore to central probing rate when app limitation "
            "happens and sender does not have enough packets to start four "
            "monitor intervals in PROBING");

namespace {
// Step size for rate change in PROBING mode.
const float kProbingStepSize = 0.02;   //Increase:converge fast Decrease:fluctuate small
// Base percentile step size for rate change in DECISION_MADE mode.
const float kDecisionMadeStepSize = 0.02f; //default:0.02
// Maximum percentile step size for rate change in DECISION_MADE mode.
const float kMaxDecisionMadeStepSize = 0.10f;
// Groups of useful monitor intervals each time in PROBING mode.
const size_t kNumIntervalGroupsInProbing = 2;
// Number of bits per byte.
const size_t kBitsPerByte = 8;
// Bandwidth filter window size in round trips.
const QuicRoundTripCount kBandwidthWindowSize = 6;
const QuicByteCount kDefaultMinimumCongestionWindow = 4 * kMaxSegmentSize;
//from PCC Vivace
const double kDelayGradientCoefficient = 900; //Vivace:b
const double kLossCoefficient = 11.35;    //Vivace:c
const double kThroughputCoefficient = 500 * 1000; 
const double kThroughputPower = 0.5;  //Vivace:t
const double kDelayGradientThreshold = 0.01;  
const double kDelayGradientNegativeBound = 10;
// The factor that converts utility gradient to sending rate change.
float kUtilityGradientToRateChangeFactor = 1.0f;
// The exponent to amplify sending rate change based on number of consecutive
// rounds in DECISION_MADE mode.
float kRateChangeAmplifyExponent = 2.0f;  //default:1.2
int incremental_rate_change_step_allowance_ = 0;
// The initial maximum rate change step size in Vivace.
float kInitialMaxStepSize = 0.05f;    //default:0.05
// The incremental rate change step size allowed on basis of current maximum
// step size every time the calculated rate change exceeds the current max.
float kIncrementalStepSize = 0.05;    //default:0.05
// The smallest sending rate change allowed by Vivace.
const QuicBandwidth kMinRateChange = QuicBandwidth::FromKBitsPerSecond(0);
// The smallest sending rate supported by Vivace.
const QuicBandwidth kMinSendingRate = QuicBandwidth::FromKBitsPerSecond(500);

const size_t kMegabit = 1000*1000;
// Number of bits per Mbit.
const size_t kKilobit = 1000;
// Number of bits per Kbit.
const size_t kampFactor = 1;//1000*1000;

}  // namespace

UtilityInfo::UtilityInfo()
    : sending_rate(QuicBandwidth::Zero()), utility(0.0) {}

UtilityInfo::UtilityInfo(QuicBandwidth rate, double utility)
    : sending_rate(rate), utility(utility) {}

PccSender::DebugState::DebugState(const PccSender& sender)
    : mode(sender.mode_),
      sending_rate(sender.interval_queue_.current().sending_rate),
      latest_rtt(sender.rtt_stats_->latest_rtt()),
      smoothed_rtt(sender.rtt_stats_->smoothed_rtt()),
      rtt_dev(sender.rtt_stats_->mean_deviation()),
      is_useful(sender.interval_queue_.current().is_useful),
      first_packet_sent_time(sender.interval_queue_.current()
                                 .first_packet_sent_time),
      last_packet_sent_time(sender.interval_queue_.current()
                                 .last_packet_sent_time),
      first_packet_number(sender.interval_queue_.current().first_packet_number),
      last_packet_number(sender.interval_queue_.current().last_packet_number),
      bytes_sent(sender.interval_queue_.current().bytes_sent),
      bytes_acked(sender.interval_queue_.current().bytes_acked),
      bytes_lost(sender.interval_queue_.current().bytes_lost),
      rtt_on_monitor_start(sender.interval_queue_.current()
                                  .rtt_on_monitor_start),
      rtt_on_monitor_end(sender.interval_queue_.current().rtt_on_monitor_end),
      latest_utility(sender.latest_utility_.utility),
      bandwidth(sender.BandwidthEstimateInner()) {}

PccSender::PccSender(const RttStats* rtt_stats,
                     const UnackedPacketMap* unacked_packets,
                     QuicPacketCount initial_congestion_window,
                     QuicPacketCount max_congestion_window,
                     Random* random,UtilityFunctionVariant fun_type)
    : mode_(STARTING),
      sending_rate_(
          QuicBandwidth::FromBytesAndTimeDelta(initial_congestion_window *
                                                   kDefaultTCPMSS,
                                               rtt_stats->initial_rtt())),
      has_seen_valid_rtt_(false),
      latest_utility_(QuicBandwidth::Zero(),0.0),
      monitor_duration_(TimeDelta::Zero()),
      direction_(INCREASE),
      rounds_(1),
      interval_queue_(rtt_stats,/*delegate=*/this),
      rtt_on_inflation_start_(TimeDelta::Zero()),
      max_cwnd_bytes_(max_congestion_window * kDefaultTCPMSS),
      min_cwnd_bytes_(kDefaultMinimumCongestionWindow),
      rtt_stats_(rtt_stats),
      unacked_packets_(unacked_packets),
      random_(random),
      max_bandwidth_(kBandwidthWindowSize, QuicBandwidth::Zero(), 0),
      last_sent_packet_(0),
      current_round_trip_end_(0),
      round_trip_count_(0),
      exit_starting_based_on_sampled_bandwidth_(
          FLAGS_exit_starting_based_on_sampled_bandwidth) {
     PccUtilityFunctionInterface *utility;
     //fun_type default is kPccUtility
     if(fun_type==kPccUtility){
        utility=new PccUtilityFunction;
     }else if(fun_type==kVivaceUtility){
        utility=new VivaceUtilityFunction(kDelayGradientCoefficient, kLossCoefficient, kThroughputCoefficient,
      kThroughputPower, kDelayGradientThreshold, kDelayGradientNegativeBound);
     }else if(fun_type==kPriorityUtility){
        utility=new PriorityUtilityFunction(kDelayGradientCoefficient, kLossCoefficient, kThroughputCoefficient,
      kThroughputPower, kDelayGradientThreshold, kDelayGradientNegativeBound);
     }else if(fun_type==kProteusUtility){
        utility=new ProteusUtilityFunction(kDelayGradientCoefficient, kLossCoefficient, kThroughputCoefficient,
      kThroughputPower, kDelayGradientThreshold, kDelayGradientNegativeBound);
     }
     utility_function_.reset(utility);
}

void PccSender::OnPacketSent(ProtoTime sent_time,
                             QuicByteCount bytes_in_flight,
                             QuicPacketNumber packet_number,
                             QuicByteCount bytes,
                             HasRetransmittableData is_retransmittable) {
  last_sent_packet_ = packet_number;

  // Do not process not retransmittable packets. Otherwise, the interval may
  // never be able to end if one of these packets gets lost.
  if (is_retransmittable != HAS_RETRANSMITTABLE_DATA) {
    return;
  }

  // Start a new monitor interval if the interval queue is empty. If latest RTT
  // is available, start a new monitor interval if (1) there is no useful
  // interval or (2) it has been more than monitor_duration since the last
  // interval starts.
  if (interval_queue_.empty() ||
      (!rtt_stats_->latest_rtt().IsZero() &&
       (interval_queue_.num_useful_intervals() == 0 ||
        sent_time - interval_queue_.current().first_packet_sent_time >
            monitor_duration_))) {
    MaybeSetSendingRate();
    // Set the monitor duration to 1.5 of min rtt.
    monitor_duration_ = TimeDelta::FromMicroseconds(
        rtt_stats_->min_rtt().ToMicroseconds() * 1.0);

    bool is_useful = CreateUsefulInterval();
    interval_queue_.EnqueueNewMonitorInterval(
        is_useful ? sending_rate_ : GetSendingRateForNonUsefulInterval(),
        is_useful, GetMaxRttFluctuationTolerance());
  }
  interval_queue_.OnPacketSent(sent_time, packet_number, bytes);
  sampler_.OnPacketSent(sent_time, packet_number, bytes, bytes_in_flight,
                          is_retransmittable);
}

void PccSender::OnCongestionEvent(bool rtt_updated,
                                  QuicByteCount bytes_in_flight,
                                  ProtoTime event_time,
                                  const AckedPacketVector& acked_packets,
                                  const LostPacketVector& lost_packets) {
  std::vector<ReceivedPacket> receipt_packets;
  UpdateBandwidthSampler(event_time, acked_packets, lost_packets,receipt_packets);
  if (!has_seen_valid_rtt_) {
    has_seen_valid_rtt_ = true;
    // Update sending rate if the actual RTT is smaller than initial rtt value
    // in RttStats, so PCC can start with larger rate and ramp up faster.
    if (rtt_stats_->latest_rtt() < rtt_stats_->initial_rtt()) {
      sending_rate_ = sending_rate_ *
          (static_cast<float>(rtt_stats_->initial_rtt().ToMicroseconds()) /
           static_cast<float>(rtt_stats_->latest_rtt().ToMicroseconds()));
    }
  }
  if (mode_ == STARTING && CheckForRttInflation()) {
    // Directly enter PROBING when rtt inflation already exceeds the tolerance
    // ratio, so as to reduce packet losses and mitigate rtt inflation.
    interval_queue_.OnRttInflationInStarting();
    EnterProbing();
    return;
  }

  interval_queue_.OnCongestionEvent(receipt_packets, lost_packets);
}

bool PccSender::CanSend(QuicByteCount bytes_in_flight) {
  if (!FLAGS_can_send_respect_congestion_window) {
    return true;
  }

  if (rtt_stats_->min_rtt() < rtt_stats_->mean_deviation()) {
    // Avoid capping bytes in flight on highly fluctuating path, because that
    // may impact throughput.
    return true;
  }
  QuicByteCount target=FLAGS_bytes_in_flight_gain * GetCongestionWindow();
  target=std::max(min_cwnd_bytes_,target);
  return bytes_in_flight <target;
}

QuicBandwidth PccSender::PacingRate(QuicByteCount bytes_in_flight) const {
  return interval_queue_.empty() ? sending_rate_
                                 : interval_queue_.current().sending_rate;
}

QuicBandwidth PccSender::BandwidthEstimateInner() const {
  return exit_starting_based_on_sampled_bandwidth_ ? max_bandwidth_.GetBest()
                                                   : QuicBandwidth::Zero();
}
QuicBandwidth PccSender::BandwidthEstimate() const {
  return interval_queue_.empty() ? sending_rate_
                                 : interval_queue_.current().sending_rate;
}
QuicByteCount PccSender::GetCongestionWindow() const {
  // Use min rtt to calculate expected congestion window except when it equals
  // 0, which happens when the connection just starts.
  return sending_rate_ * (rtt_stats_->min_rtt().IsZero()
                              ? rtt_stats_->initial_rtt()
                              : rtt_stats_->min_rtt());
}

bool PccSender::InSlowStart() const { return false; }

bool PccSender::InRecovery() const { return false; }

bool PccSender::ShouldSendProbingPacket() const { return false; }

QuicByteCount PccSender::GetSlowStartThreshold() const { return 0; }

CongestionControlType PccSender::GetCongestionControlType() const {
  return kPCC;
}

void PccSender::OnApplicationLimited(QuicByteCount bytes_in_flight) {
  if (!exit_starting_based_on_sampled_bandwidth_ ||
      bytes_in_flight >= GetCongestionWindow()) {
    return;
  }
  sampler_.OnAppLimited();
}

std::string PccSender::GetDebugState() const {
  if (interval_queue_.empty()) {
    return "pcc??";
  }

  std::ostringstream stream;
  stream << ExportDebugState();
  return stream.str();
}

PccSender::DebugState PccSender::ExportDebugState() const {
  return DebugState(*this);
}
QuicBandwidth PccSender::GetMinRate(){
    TimeDelta rtt=(rtt_stats_->smoothed_rtt().IsZero()
                              ? rtt_stats_->initial_rtt()
                              : rtt_stats_->smoothed_rtt());
    return QuicBandwidth::FromBytesAndTimeDelta(min_cwnd_bytes_,rtt);
}
void PccSender::UpdateBandwidthSampler(ProtoTime event_time,
                                       const AckedPacketVector& acked_packets,
                                       const LostPacketVector& lost_packets,
                                       ReceivedPacketVector&receipt_packets) {
  // This function should not be called if latched value of
  // FLAGS_exit_starting_based_on_sampled_bandwidth is false.
  DCHECK(exit_starting_based_on_sampled_bandwidth_);

  // Update round trip count if largest acked packet number exceeds largest
  // packet number in current round trip.
  if (!acked_packets.empty() &&
      acked_packets.rbegin()->packet_number > current_round_trip_end_) {
    round_trip_count_++;
    current_round_trip_end_ = last_sent_packet_;
  }
  // Calculate bandwidth based on the acked packets.
  for (const AckedPacket& packet : acked_packets) {
    if (packet.bytes_acked == 0) {
      receipt_packets.emplace_back(packet.packet_number,packet.bytes_acked);
      continue;
    }
    BandwidthSample bandwidth_sample =
        sampler_.OnPacketAcknowledged(event_time, packet.packet_number);
    receipt_packets.emplace_back(bandwidth_sample.sent_time,packet.packet_number,
    bandwidth_sample.rtt,packet.bytes_acked);
    if (!bandwidth_sample.state_at_send.is_app_limited||
        bandwidth_sample.bandwidth > BandwidthEstimateInner()) {
      max_bandwidth_.Update(bandwidth_sample.bandwidth, round_trip_count_);
    }
  }
  // Remove lost and obsolete packets from bandwidth sampler.
  for (const LostPacket& packet : lost_packets) {
    sampler_.OnPacketLost(packet.packet_number);
  }
  sampler_.RemoveObsoletePackets(unacked_packets_->GetLeastUnacked());
}

void PccSender::OnUtilityAvailable(
    const std::vector<const MonitorInterval *>& useful_intervals) {
  // Calculate the utilities for all available intervals.
  std::vector<UtilityInfo> utility_info;
  std::cout << std::endl;
  for(size_t i = 0; i < useful_intervals.size(); ++i) {
    utility_info.push_back(
        UtilityInfo(useful_intervals[i]->sending_rate,
                    utility_function_->CalculateUtility(useful_intervals[i])));
  }
  std::cout << "state now:" << PccSenderModeToString(mode_) << std::endl;
  // const QuicBandwidth kLargeSendingRate = QuicBandwidth::FromKBitsPerSecond(2700);
  // if utility < 0, set sending rate to 0.90*sending rate.
  // if (utility_info[0].utility < 0 && sending_rate_ > kLargeSendingRate ) {
  //   sending_rate_ = sending_rate_ * 0.90;
  //   mode_ = PROBING;
  //   std::cout << "Utility < 0,decrease sengding rate,change sengding rate to:" << sending_rate_ << std::endl;
  //   std::cout << "And change mode to" << PccSenderModeToString(mode_) << std::endl;
  //   return;
  // }

  switch (mode_) {
    case STARTING:
      DCHECK_EQ(1u, utility_info.size());
      if (utility_info[0].utility > latest_utility_.utility) {
        // Stay in STARTING mode. Double the sending rate and update
        // latest_utility.
        sending_rate_ = sending_rate_ * 2;
        latest_utility_ = utility_info[0];
        ++rounds_;
      } else {
        // Enter PROBING mode if utility decreases.
        std::cout << "In STARTING Enter Probing." << std::endl;
        EnterProbing();
      }
      break;
    case PROBING:
      if (CanMakeDecision(utility_info)) {
        if (FLAGS_restore_central_rate_upon_app_limited &&
            interval_queue_.current().is_useful) {
          // If there is no non-useful interval in this round of PROBING, sender
          // needs to change sending_rate_ back to central rate.
          RestoreCentralSendingRate();
        }
        DCHECK_EQ(2 * kNumIntervalGroupsInProbing, utility_info.size());
        // Enter DECISION_MADE mode if a decision is made.
        direction_ = (utility_info[0].utility > utility_info[1].utility)
                         ? ((utility_info[0].sending_rate >
                             utility_info[1].sending_rate)
                                ? INCREASE
                                : DECREASE)
                         : ((utility_info[0].sending_rate >
                             utility_info[1].sending_rate)
                                ? DECREASE
                                : INCREASE);
        // latest_utility_ =
        //     std::max(utility_info[2 * kNumIntervalGroupsInProbing - 2].,
        //              utility_info[2 * kNumIntervalGroupsInProbing - 1].utility);
        std::cout << "In Probing Enter DecisionMade." << std::endl;
        EnterDecisionMade(utility_info);
      } else {
        // Stays in PROBING mode.
        std::cout << "In Probing Enter Probing." << std::endl;
        EnterProbing();
      }
      if ((rounds_ > 1 || mode_ == DECISION_MADE) && sending_rate_ <= kMinSendingRate) {
        sending_rate_ = kMinSendingRate;
        incremental_rate_change_step_allowance_ = 0;
        rounds_ = 1;
        mode_ = STARTING;
      }
      break;
    case DECISION_MADE:
      DCHECK_EQ(1u, utility_info.size());
      ////PCC : Decision_Made
      // if (utility_info[0].utility > latest_utility_) {
      //   // Remain in DECISION_MADE mode. Keep increasing or decreasing the
      //   // sending rate.
      //   ++rounds_;
      //   if (direction_ == INCREASE) {
      //     sending_rate_ = sending_rate_ *
      //                     (1 + std::min(rounds_ * kDecisionMadeStepSize,
      //                                   kMaxDecisionMadeStepSize));
      //   } else {
      //     sending_rate_ = sending_rate_ *
      //                     (1 - std::min(rounds_ * kDecisionMadeStepSize,
      //                                   kMaxDecisionMadeStepSize));
      //   }
      //   latest_utility_ = utility_info[0];
      if ((direction_ == INCREASE &&
           utility_info[0].utility > latest_utility_.utility &&
           utility_info[0].sending_rate > latest_utility_.sending_rate) ||
          (direction_ == INCREASE &&
           utility_info[0].utility < latest_utility_.utility &&
           utility_info[0].sending_rate < latest_utility_.sending_rate) ||
          (direction_ == DECREASE &&
           utility_info[0].utility > latest_utility_.utility &&
           utility_info[0].sending_rate < latest_utility_.sending_rate) ||
          (direction_ == DECREASE &&
           utility_info[0].utility < latest_utility_.utility &&
           utility_info[0].sending_rate > latest_utility_.sending_rate)) {
        // Remain in DECISION_MADE mode. Keep increasing or decreasing the
        // sending rate.
        std::cout << "In DecisionMade Enter DecisionMade." << std::endl;
        EnterDecisionMade(utility_info);
        latest_utility_ = utility_info[0];
      } else {
        // Enter PROBING mode if utility decreases.
        std::cout << "In DecisionMade Enter Probing_Vivace." << std::endl;
        std::cout << "difference happen between direction and utility." << std::endl;
        EnterProbing_vivace(utility_info);
      }
      break;
  }
}

bool PccSender::CreateUsefulInterval() const {
  if (rtt_stats_->smoothed_rtt().ToMicroseconds() == 0) {
    // Create non useful intervals upon starting a connection, until there is
    // valid rtt stats.
    QUIC_BUG_IF(mode_ != STARTING);
    return false;
  }
  // In STARTING and DECISION_MADE mode, there should be at most one useful
  // intervals in the queue; while in PROBING mode, there should be at most
  // 2 * kNumIntervalGroupsInProbing.
  size_t max_num_useful =
      (mode_ == PROBING) ? 2 * kNumIntervalGroupsInProbing : 1;
  return interval_queue_.num_useful_intervals() < max_num_useful;
}

QuicBandwidth PccSender::GetSendingRateForNonUsefulInterval() const {
  switch (mode_) {
    case STARTING:
      // Use halved sending rate for non-useful intervals in STARTING.
      return sending_rate_ * 0.5;
    case PROBING:
      // Use the smaller probing rate in PROBING.
      return sending_rate_ * (1 - kProbingStepSize);
    case DECISION_MADE:
      // Use the last (smaller) sending rate if the sender is increasing sending
      // rate in DECISION_MADE. Otherwise, use the current sending rate.
      return direction_ == DECREASE
          ? sending_rate_
          : sending_rate_ *
                (1.0 / (1 + std::min(rounds_ * kDecisionMadeStepSize,
                                     kMaxDecisionMadeStepSize)));
  }
}

void PccSender::MaybeSetSendingRate() {
  if (mode_ != PROBING || (interval_queue_.num_useful_intervals() ==
                               2 * kNumIntervalGroupsInProbing &&
                           !interval_queue_.current().is_useful)) {
    // Do not change sending rate when (1) current mode is STARTING or
    // DECISION_MADE (since sending rate is already changed in
    // OnUtilityAvailable), or (2) more than 2 * kNumIntervalGroupsInProbing
    // intervals have been created in PROBING mode.
    return;
  }

  if (interval_queue_.num_useful_intervals() != 0) {
    // Restore central sending rate.
    std::cout << "From MaybeSetSendingRate to RestorCentralSendingRate." << std::endl;
    RestoreCentralSendingRate();

    if (interval_queue_.num_useful_intervals() ==
        2 * kNumIntervalGroupsInProbing) {
      // This is the first not useful monitor interval, its sending rate is the
      // central rate.
      return;
    }
  }

  // Sender creates several groups of monitor intervals. Each group comprises an
  // interval with increased sending rate and an interval with decreased sending
  // rate. Which interval goes first is randomly decided.
  if (interval_queue_.num_useful_intervals() % 2 == 0) {
    direction_ = (random_->nextInt() % 2 == 1) ? INCREASE : DECREASE;
  } else {
    direction_ = (direction_ == INCREASE) ? DECREASE : INCREASE;
  }
  if (direction_ == INCREASE) {
    sending_rate_ = sending_rate_ * (1 + kProbingStepSize);
  } else {
    sending_rate_ = sending_rate_ * (1 - kProbingStepSize);
  }
  std::cout << "In MayBeSetSendingRate phase,change sending rate to: " << sending_rate_ << std::endl << std::endl;
  // sending_rate_=std::max(sending_rate_,GetMinRate());
}

float PccSender::GetMaxRttFluctuationTolerance() const {
  if (mode_ == PROBING) {
    // No rtt fluctuation tolerance during PROBING.
    return 0.0f;
  }

  float tolerance_ratio =
      mode_ == STARTING
          ? FLAGS_max_rtt_fluctuation_tolerance_ratio_in_starting
          : FLAGS_max_rtt_fluctuation_tolerance_ratio_in_decision_made;

  if (FLAGS_enable_rtt_deviation_based_early_termination) {
    float tolerance_gain =
        mode_ == STARTING
            ? FLAGS_rtt_fluctuation_tolerance_gain_in_starting
            : FLAGS_rtt_fluctuation_tolerance_gain_in_decision_made;
    tolerance_ratio = std::min(
        tolerance_ratio,
        tolerance_gain *
            static_cast<float>(rtt_stats_->mean_deviation().ToMicroseconds()) /
            static_cast<float>(rtt_stats_->SmoothedOrInitialRtt()
                                   .ToMicroseconds()));
  }

  return tolerance_ratio;
}

bool PccSender::CanMakeDecision(
    const std::vector<UtilityInfo>& utility_info) const {
  // Determine whether increased or decreased probing rate has better utility.
  // Cannot make decision if number of utilities are less than
  // 2 * kNumIntervalGroupsInProbing. This happens when sender does not have
  // enough data to send.
  if (utility_info.size() < 2 * kNumIntervalGroupsInProbing) {
    return false;
  }

  bool increase = false;
  // All the probing groups should have consistent decision. If not, directly
  // return false.
  for (size_t i = 0; i < kNumIntervalGroupsInProbing; ++i) {
    bool increase_i =
        utility_info[2 * i].utility > utility_info[2 * i + 1].utility
            ? utility_info[2 * i].sending_rate >
                  utility_info[2 * i + 1].sending_rate
            : utility_info[2 * i].sending_rate <
                  utility_info[2 * i + 1].sending_rate;

    if (i == 0) {
      increase = increase_i;    //first interval vs other interval
    }
    // Cannot make decision if groups have inconsistent results.
    if (increase_i != increase) {
      return false;
    }
  }

  return true;
}

void PccSender::EnterProbing() {
  switch (mode_) {
    case STARTING:
      // Fall back to the minimum between halved sending rate and
      // max bandwidth * (1 - 0.05) if there is valid bandwidth sample.
      // Otherwise, simply halve the current sending rate.
      sending_rate_ = sending_rate_ * 0.5;
      if (!BandwidthEstimateInner().IsZero()) {
        DCHECK(exit_starting_based_on_sampled_bandwidth_);
        sending_rate_ = std::min(sending_rate_,
                                 BandwidthEstimateInner() * (1 - kProbingStepSize));
      }
      break;
    case DECISION_MADE:
      FALLTHROUGH_INTENDED;
    case PROBING:
      // Reset sending rate to central rate when sender does not have enough
      // data to send more than 2 * kNumIntervalGroupsInProbing intervals.
      std::cout << "From Probing to Probing." << std::endl;
      RestoreCentralSendingRate();
      break;
  }

  if (mode_ == PROBING) {
    ++rounds_;
    return;
  }

  mode_ = PROBING;
  std::cout << "chane state to:" << PccSenderModeToString(mode_) << std::endl;
  rounds_ = 1;
}

void PccSender::EnterProbing_vivace (const std::vector<UtilityInfo>& utility_info) {
  assert(DECISION_MADE == mode_);

  rounds_ = 1;

  std::cout << "Probing vivace to compute rate change." << std::endl; 
  QuicBandwidth rate_change = ComputeRateChange(utility_info);
  std::cout << "rate change:" << rate_change << std::endl;
  if (direction_ == INCREASE) {
    sending_rate_ = sending_rate_ - rate_change;
  } else {
    sending_rate_ = sending_rate_ +  rate_change;
  }
  std::cout << "direction:" << DirectionToString(direction_) << std::endl;
  std::cout << "sending rate:" << sending_rate_ << std::endl << std::endl;

  if (sending_rate_ < kMinSendingRate) {
    sending_rate_ = kMinSendingRate;
    incremental_rate_change_step_allowance_ = 0;
  }
  mode_ = PROBING;
  std::cout << "chane state to:" << PccSenderModeToString(mode_) << std::endl;
  

}

void PccSender::EnterDecisionMade(const std::vector<UtilityInfo>& utility_info) {
  DCHECK_EQ(PROBING, mode_);

  // Change sending rate from central rate based on the probing rate with higher
  // utility.
  // if (direction_ == INCREASE) {
  //   sending_rate_ = sending_rate_ * (1 + kProbingStepSize) *
  //                   (1 + kDecisionMadeStepSize);
  // } else {
  //   sending_rate_ = sending_rate_ * (1 - kProbingStepSize) *
  //                   (1 - kDecisionMadeStepSize);
  // }
  // sending_rate_=std::max(sending_rate_,GetMinRate());
  // mode_ = DECISION_MADE;
  // rounds_ = 1;

  //  PCC Vivace Decision Made
  if (mode_ == PROBING) {
    if (direction_ == INCREASE) {
      sending_rate_ = sending_rate_ * (1 + kProbingStepSize);
    } else{
      sending_rate_ = sending_rate_ * (1 - kProbingStepSize);
    }
  }

  rounds_ = mode_ == PROBING ? 1 : rounds_ + 1;

  QuicBandwidth rate_change = ComputeRateChange(utility_info);
  std::cout << "rate change:" << rate_change << std::endl;
  if (direction_ == INCREASE) {
    sending_rate_ = sending_rate_ + rate_change;
  } else {
    sending_rate_ = sending_rate_ -  rate_change;
  }
  std::cout << "direction:" << DirectionToString(direction_) << std::endl;
  std::cout << "sending rate:" << sending_rate_ << std::endl << std::endl;

  if (sending_rate_ < GetMinRate()) {
    sending_rate_ = GetMinRate();
    std::cout << "sending rate = min rate:" << GetMinRate() << std::endl;
    mode_ = PROBING;
    std::cout << "chane state to:" << PccSenderModeToString(mode_) << std::endl;
    rounds_ = 1;
    incremental_rate_change_step_allowance_ = 0;
  } else {
    mode_ = DECISION_MADE;
    std::cout << "chane state to:" << PccSenderModeToString(mode_) << std::endl;
  }
}

void PccSender::RestoreCentralSendingRate() {
  switch (mode_) {
    case STARTING:
      // The sending rate upon exiting STARTING is set separately. This function
      // should not be called while sender is in STARTING mode.
      QUIC_BUG << "Attempt to set probing rate while in STARTING";
      break;
    case PROBING:
      // Change sending rate back to central probing rate.
      if (interval_queue_.current().is_useful) {
        std::cout << "In RestoreCentralSendingRate, change rate."<<std::endl;
        if (direction_ == INCREASE) {
          sending_rate_ = sending_rate_ * (1.0 / (1 + kProbingStepSize));
        } else {
          sending_rate_ = sending_rate_ * (1.0 / (1 - kProbingStepSize));
        }
      }
      std::cout << "In Probing RestoreCentralSendingRate,bake rate to:" << sending_rate_ << std::endl << std::endl;
      break;
    case DECISION_MADE:
      if (direction_ == INCREASE) {
        sending_rate_ = sending_rate_ *
                        (1.0 / (1 + std::min(rounds_ * kDecisionMadeStepSize,
                                             kMaxDecisionMadeStepSize)));
      } else {
        sending_rate_ = sending_rate_ *
                        (1.0 / (1 - std::min(rounds_ * kDecisionMadeStepSize,
                                             kMaxDecisionMadeStepSize)));
      }
      std::cout << "Decision_Made Enter RestoreCentralSendingRate,bake rate to:" << sending_rate_ << std::endl << std::endl;
      break;
  }
}

bool PccSender::CheckForRttInflation() {
  if (interval_queue_.empty() ||
      interval_queue_.front().rtt_on_monitor_start.IsZero() ||
      rtt_stats_->latest_rtt() <= rtt_stats_->smoothed_rtt()) {
    // RTT is not inflated if latest RTT is no larger than smoothed RTT.
    rtt_on_inflation_start_ = TimeDelta::Zero();
    return false;
  }

  // Once the latest RTT exceeds the smoothed RTT, store the corresponding
  // smoothed RTT as the RTT at the start of inflation. RTT inflation will
  // continue as long as latest RTT keeps being larger than smoothed RTT.
  if (rtt_on_inflation_start_.IsZero()) {
    rtt_on_inflation_start_ = rtt_stats_->smoothed_rtt();
  }

  const float max_inflation_ratio = 1 + GetMaxRttFluctuationTolerance();
  const TimeDelta rtt_on_monitor_start =
      FLAGS_trigger_early_termination_based_on_interval_queue_front
          ? interval_queue_.front().rtt_on_monitor_start
          : interval_queue_.current().rtt_on_monitor_start;
  bool is_inflated =
      max_inflation_ratio * rtt_on_monitor_start < rtt_stats_->smoothed_rtt();
  if (!is_inflated &&
      FLAGS_enable_early_termination_based_on_latest_rtt_trend) {
    // If enabled, check if early termination should be triggered according to
    // the stored smoothed rtt on inflation start.
    is_inflated = max_inflation_ratio * rtt_on_inflation_start_ <
                      rtt_stats_->smoothed_rtt();
  }
  if (is_inflated) {
    // RTT is inflated by more than the tolerance, and early termination will be
    // triggered. Reset the rtt on inflation start.
    rtt_on_inflation_start_ = TimeDelta::Zero();
  }
  return is_inflated;
}

QuicBandwidth PccSender::ComputeRateChange(
  const std::vector<UtilityInfo>& utility_info) const {
  assert(mode_ != STARTING);
  
  //Compute rate difference between higher and lower sending rate,
  //as well as their utility differemce
  QuicBandwidth delta_sending_rate = QuicBandwidth::Zero();
  float delta_utility = 0.0;
  if (mode_ == PROBING) {
    delta_sending_rate = 
        std::max(utility_info[0].sending_rate, utility_info[1].sending_rate) -
        std::min(utility_info[0].sending_rate, utility_info[1].sending_rate);

    for (size_t i = 0; i < kNumIntervalGroupsInProbing; ++i) {
      bool increase_i =
          utility_info[2 * i].utility > utility_info[2 * i + 1].utility
              ? utility_info[2 * i].sending_rate >
                    utility_info[2 * i + 1].sending_rate
              : utility_info[2 * i].sending_rate <
                    utility_info[2 * i + 1].sending_rate;
      if ((increase_i && direction_ == DECREASE) || 
          (!increase_i && direction_ == INCREASE)){
            continue;
          }
      delta_utility = delta_utility +
          std::max(utility_info[2 * i].utility,utility_info[2 * i +1].utility) -
          std::min(utility_info[2 * i].utility,utility_info[2 * i +1].utility);
  }
  delta_utility /= kNumIntervalGroupsInProbing;
}else{
    delta_sending_rate = 
      std::max(utility_info[0].sending_rate, latest_utility_.sending_rate) -
      std::min(utility_info[0].sending_rate, latest_utility_.sending_rate);

    delta_utility = 
      std::max(utility_info[0].utility, latest_utility_.utility) -
      std::min(utility_info[0].utility, latest_utility_.utility);   //diff
}

  assert(!delta_sending_rate.IsZero());

  float utility_gradient =
      delta_utility / delta_sending_rate.ToBitsPerSecond();  //ToBitPerSecond?
  
  std::cout << utility_gradient << std::endl;

  QuicBandwidth rate_change = QuicBandwidth::FromBitsPerSecond(
      utility_gradient * kampFactor * kMegabit * kUtilityGradientToRateChangeFactor);  //kMegabit
  std::cout << "compute1 rate change:" << rate_change << std::endl;
    
//PCC QUIC——meng tong
  if (mode_ == DECISION_MADE) {
    // Amplify rate change amount when sending rate changes towards the same
    // direction more than once.
    rate_change = rate_change * pow(static_cast<float>((rounds_ + 1) / 2),
                                    kRateChangeAmplifyExponent);
  } else {
    // Reset allowed incremental rate change step size upon entering PROBING.
    incremental_rate_change_step_allowance_ = 0;
  }

  QuicBandwidth max_allowed_rate_change =
      sending_rate_ * (kInitialMaxStepSize +
                       kIncrementalStepSize /2 * static_cast<float>(
                           incremental_rate_change_step_allowance_));
  if (rate_change > max_allowed_rate_change) {
    rate_change = max_allowed_rate_change;
    // Increase incremental rate change step size if the calculated rate change
    // exceeds the current maximum.
    incremental_rate_change_step_allowance_++;
  } else if (incremental_rate_change_step_allowance_ > 0) {
    // Reduce incremental rate change allowance if calculated rate is smaller
    // than the current maximum.
    incremental_rate_change_step_allowance_--;
  }
  std::cout << "compute rate change:" << rate_change << std::endl;
  return std::max(rate_change ,kMinRateChange);
}

static std::string PccSenderModeToString(PccSender::SenderMode mode) {
  switch (mode) {
    case PccSender::STARTING:
      return "STARTING";
    case PccSender::PROBING:
      return "PROBING";
    case PccSender::DECISION_MADE:
      return "DECISION_MADE";
  }
  return "???";
}

static std::string DirectionToString(PccSender::RateChangeDirection direction) {
  switch (direction) {
    case PccSender::INCREASE:
      return "INCREASE";
    case PccSender::DECREASE:
      return "DECREASE";
  }
  return "unknow!!!!";
}

std::ostream& operator<<(std::ostream& os, const PccSender::DebugState& state) {
  os << "Mode: " << PccSenderModeToString(state.mode) << std::endl;
  os << "Sending rate: " << state.sending_rate.ToKBitsPerSecond() << std::endl;
  os << "Latest rtt: " << state.latest_rtt.ToMicroseconds() << std::endl;
  os << "Smoothed rtt: " << state.smoothed_rtt.ToMicroseconds() << std::endl;
  os << "Rtt deviation: " << state.rtt_dev.ToMicroseconds() << std::endl;
  os << "Monitor useful: " << (state.is_useful ? "yes" : "no") << std::endl;
  os << "Monitor packet sent time: "
     << state.first_packet_sent_time.ToDebuggingValue() << " -> "
     << state.last_packet_sent_time.ToDebuggingValue() << std::endl;
  os << "Monitor packet number: " << state.first_packet_number << " -> "
     << state.last_packet_number << std::endl;
  os << "Monitor bytes: " << state.bytes_sent << " (sent), "
     << state.bytes_acked << " (acked), " << state.bytes_lost << " (lost)"
     << std::endl;
  os << "Monitor rtt change: " << state.rtt_on_monitor_start.ToMicroseconds()
     << " -> " << state.rtt_on_monitor_end.ToMicroseconds() << std::endl;
  os << "Latest utility: " << state.latest_utility << std::endl;
  os << "Bandwidth sample: " << state.bandwidth.ToKBitsPerSecond() << std::endl;

  return os;
}

}  // namespace qcc
