#pragma once

#include "Action.h"

class MicromorphicContinuumAction : public Action
{
public:
  static InputParameters validParams();

  MicromorphicContinuumAction( const InputParameters & params );

  void act();

protected:
  void addKernels();
  void addMaterial();

  /// these parameters are not passed to the invoked kernels
  const static std::vector< std::string > excludedParameters;

  const unsigned int _ndisp;
  const unsigned int _nmicro_disp_grad;
};
