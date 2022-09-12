//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COMPUTEPROJECTOR_H
#define COMPUTEPROJECTOR_H

#include "AuxKernel.h"
#include "NodalOverlapUserObject.h"
#include "ElementIntegrateUserObject.h"
#include "occonfiguration.h"
#include "overlap_coupling.h"

// Forward Declarations
class ComputeProjectorAux;

template <>
InputParameters validParams<ComputeProjectorAux>();

/**
 * Coupled auxiliary value
 */
class ComputeProjectorAux : public AuxKernel
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  ComputeProjectorAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const NodalOverlapUserObject & _nodal_overlap;
  const ElementIntegrateUserObject & _element_integrate;

};

#endif // COMPUTEPROJECTOR_H
