//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef DUMMYAUX_H
#define DUMMYAUX_H

#include "AuxKernel.h"

// Forward Declarations
class DummyAux;

template <>
InputParameters validParams<DummyAux>();

/**
 * Coupled auxiliary value
 */
class DummyAux : public AuxKernel
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  DummyAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const VariableValue & _priorvar;

};

#endif // DUMMYAUX_H
