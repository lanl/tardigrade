//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ELEMENTUOAUX_H
#define ELEMENTUOAUX_H

#include "AuxKernel.h"
#include "ElementIntegrateUserObject.h"

// Forward Declarations
class ElementUOAux;

template <>
InputParameters validParams<ElementUOAux>();

/**
 * Coupled auxiliary value
 */
class ElementUOAux : public AuxKernel
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  ElementUOAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const VariableValue & _priorvar;

  const ElementIntegrateUserObject & _element_integrate;

};

#endif // ELEMENTUOAUX_H
