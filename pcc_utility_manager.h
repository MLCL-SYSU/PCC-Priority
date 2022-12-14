#pragma once
#include "pcc_monitor_interval_queue.h"
namespace qcc {
class PccUtilityFunctionInterface {
public:
  virtual double CalculateUtility(const MonitorInterval* interval) const = 0;
  virtual ~PccUtilityFunctionInterface() {};
};
class PccUtilityFunction:public PccUtilityFunctionInterface{
public:
    PccUtilityFunction(){}
    ~PccUtilityFunction() override{}
    double CalculateUtility(const MonitorInterval* interval) const override;
};
class VivaceUtilityFunction:public PccUtilityFunctionInterface{
public:
    VivaceUtilityFunction(double delay_gradient_coefficient,
                        double loss_coefficient,
                        double throughput_coefficient,
                        double throughput_power,
                        double delay_gradient_threshold,
                        double delay_gradient_negative_bound);
    ~VivaceUtilityFunction() override{}
    double CalculateUtility(const MonitorInterval* monitor_interval) const override;
private:
  const double delay_gradient_coefficient_;
  const double loss_coefficient_;
  const double throughput_power_;
  const double throughput_coefficient_;
  const double delay_gradient_threshold_;
  const double delay_gradient_negative_bound_; 
};

class PriorityUtilityFunction:public PccUtilityFunctionInterface{
public:
    PriorityUtilityFunction(double delay_gradient_coefficient,
                        double loss_coefficient,
                        double throughput_coefficient,
                        double throughput_power,
                        double delay_gradient_threshold,
                        double delay_gradient_negative_bound);
    ~PriorityUtilityFunction() override{}
    double CalculateUtility(const MonitorInterval* monitor_interval) const override;
private:
  const double delay_gradient_coefficient_;
  const double loss_coefficient_;
  const double throughput_power_;
  const double throughput_coefficient_;
  const double delay_gradient_threshold_;
  const double delay_gradient_negative_bound_; 
};

class ProteusUtilityFunction:public PccUtilityFunctionInterface{
public:
    ProteusUtilityFunction(double delay_gradient_coefficient,
                        double loss_coefficient,
                        double throughput_coefficient,
                        double throughput_power,
                        double delay_gradient_threshold,
                        double delay_gradient_negative_bound);
    ~ProteusUtilityFunction() override{}
    double CalculateUtility(const MonitorInterval* monitor_interval) const override;
private:
  const double delay_gradient_coefficient_;
  const double loss_coefficient_;
  const double throughput_power_;
  const double throughput_coefficient_;
  const double delay_gradient_threshold_;
  const double delay_gradient_negative_bound_; 
};
}  // namespace qcc
