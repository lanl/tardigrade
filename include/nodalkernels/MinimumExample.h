#pragma once

#include "NodalKernel.h"
#include "OverlapCoupling.h"

// Forward Declarations
class MinimumExample;

template <>
InputParameters validParams<MinimumExample>();

/**
 * Represents the rate in a simple ODE of du/dt = rate
 */
class MinimumExample : public NodalKernel
{
public:
  /**
   * Constructor initializes the rate
   */
  static InputParameters validParams();

  MinimumExample(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

};
