#include "MinimumExample.h"

registerMooseObject("MooseApp", MinimumExample);

InputParameters
MinimumExample::validParams()
{
  InputParameters params = NodalKernel::validParams();
  params.addClassDescription("Minimum example of unexpected nodal kernel behavior");
  return params;
}

MinimumExample::MinimumExample(const InputParameters & parameters)
  : NodalKernel(parameters)
{
}

Real
MinimumExample::computeQpResidual()
{
  return (*_current_node)(2);
}

Real
MinimumExample::computeQpJacobian()
{
  return 0;
}
