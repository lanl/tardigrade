//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ProjectionUserObject.h"

#include "libmesh/quadrature.h"

registerMooseObject("tardigradeApp", ProjectionUserObject);

template <>
InputParameters
validParams<ProjectionUserObject>()
{
    InputParameters params = validParams<ElementUserObject>();
    params.addRequiredParam<SubdomainName>("macroscale_domain", "The micromorphic macroscale domain");
    return params;
}

ProjectionUserObject::ProjectionUserObject(const InputParameters & parameters)
    : ElementUserObject(parameters),
    _macroscale_domain(getParam<SubdomainName>("macroscale_domain"))
  
{
}

void
ProjectionUserObject::initialize()
{
    _console << "Initializing Projection UserObject " << name() << std::endl;
}

void
ProjectionUserObject::execute()
{
  
}

void
ProjectionUserObject::threadJoin(const UserObject & y)
{
    const ProjectionUserObject & pps = static_cast<const ProjectionUserObject &>(y);
}
