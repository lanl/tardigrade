//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NodalUOAux.h"

registerMooseObject("tardigradeApp", NodalUOAux);

template <>
InputParameters
validParams<NodalUOAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addCoupledVar("priorvar", 0, "A variable which should be computed prior to the current AuxKernel.");
//  params.addRequiredParam<UserObjectName>("nodal_overlap_userobject", "The NodalOverlapUserObject that detects which micro-scale nodes overlap with the macro-scale");
  params.addRequiredParam<UserObjectName>("nodal_userobject", "The NodalUserObject to force execution");

  return params;
}

NodalUOAux::NodalUOAux(const InputParameters & parameters)
  : AuxKernel(parameters),
  _priorvar(coupledValue("priorvar")),
//  _nodal_overlap(getUserObject<NodalOverlapUserObject>("nodal_overlap_userobject"))
  _nodal_userobject(getUserObject<NodalUserObject>("nodal_userobject"))
{
}

/**
 * Auxiliary Kernels override computeValue() instead of computeQpResidual().  Aux Variables
 * are calculated either one per elemenet or one per node depending on whether we declare
 * them as "Elemental (Constant Monomial)" or "Nodal (First Lagrange)".  No changes to the
 * source are necessary to switch from one type or the other.
 */

Real
NodalUOAux::computeValue()
{
    return 0;
}
