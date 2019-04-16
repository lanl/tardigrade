//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef NODALUOAUX_H
#define NODALUOAUX_H

#include "AuxKernel.h"
//#include "NodalOverlapUserObject.h"
#include "NodalUserObject.h"

// Forward Declarations
class NodalUOAux;

template <>
InputParameters validParams<NodalUOAux>();

/**
 * Coupled auxiliary value
 */
class NodalUOAux : public AuxKernel
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  NodalUOAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const VariableValue & _priorvar;

  const NodalUserObject & _nodal_userobject;
//  const NodalOverlapUserObject & _nodal_overlap;

};

#endif // NODALUOAUX_H
