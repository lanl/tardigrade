//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementOverlapUserObject.h"

#include "libmesh/bounding_box.h"
#include "libmesh/mesh_tools.h"

registerMooseObject("tardigradeApp", ElementOverlapUserObject);

template <>
InputParameters
validParams<ElementOverlapUserObject>()
{
    InputParameters params = validParams<ElementUserObject>();
    return params;
}

ElementOverlapUserObject::ElementOverlapUserObject(const InputParameters & parameters)
    : ElementUserObject(parameters)
{
}

void
ElementOverlapUserObject::initialize()
{
    _console << "Initializing Element Overlap UserObject: " << name() << std::endl;
}

void
ElementOverlapUserObject::execute()
{

}

void
ElementOverlapUserObject::threadJoin(const UserObject & y)
{

}

void
ElementOverlapUserObject::finalize()
{
    std::cout << "Finalizing Element Overlap\n";
}
