//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementUOAux.h"

registerMooseObject("tardigradeApp", ElementUOAux);

template <>
InputParameters
validParams<ElementUOAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addCoupledVar("priorvar", 0, "A variable which should be computed prior to the current AuxKernel.");
  params.addRequiredParam<UserObjectName>("element_integrate_userobject", "The ElementIntegrateUserObject that computes shape-function weighted integrals at the nodes");
  return params;
}

ElementUOAux::ElementUOAux(const InputParameters & parameters)
  : AuxKernel(parameters),
  _priorvar(coupledValue("priorvar")),
  _element_integrate(getUserObject<ElementIntegrateUserObject>("element_integrate_userobject"))
{
}

/**
 * Auxiliary Kernels override computeValue() instead of computeQpResidual().  Aux Variables
 * are calculated either one per elemenet or one per node depending on whether we declare
 * them as "Elemental (Constant Monomial)" or "Element (First Lagrange)".  No changes to the
 * source are necessary to switch from one type or the other.
 */

Real
ElementUOAux::computeValue()
{
    return 0;
}
