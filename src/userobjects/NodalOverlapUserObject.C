//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NodalOverlapUserObject.h"

#include "libmesh/bounding_box.h"
#include "libmesh/mesh_tools.h"

registerMooseObject("tardigradeApp", NodalOverlapUserObject);

template <>
InputParameters
validParams<NodalOverlapUserObject>()
{
    InputParameters params = validParams<NodalUserObject>();
    params.addRequiredParam<SubdomainName>("macroscale_domain", "The micromorphic macroscale domain");
    return params;
}

NodalOverlapUserObject::NodalOverlapUserObject(const InputParameters & parameters)
    : NodalUserObject(parameters),
    _macroscale_domain(getParam<SubdomainName>("macroscale_domain"))
  
{
}

void
NodalOverlapUserObject::initialize()
{
    _console << "Initializing Nodal Overlap UserObject " << name() << std::endl;

    //!Set the id of the macro-scale domain
    _macro_id = _mesh.getSubdomainID(_macroscale_domain);
    
    //!Create a bounding box that encompasses the macro-scale domain
    _macro_bounding_box = MeshTools::create_subdomain_bounding_box(_mesh, _macro_id);

    //!Output the bounds to the terminal
    std::cout << "Macro-scale bounding box:\n";
    std::cout << "minimum: ";
    _macro_bounding_box.min().write_unformatted(std::cout);
    std::cout << "maximum: ";
    _macro_bounding_box.max().write_unformatted(std::cout);
}

void
NodalOverlapUserObject::execute()
{
    //!Check if the node is located within the macro-scale's bounding box
    if (_macro_bounding_box.contains_point(_current_node))
    {
        
    }
}

void
NodalOverlapUserObject::threadJoin(const UserObject & y)
{
    const NodalOverlapUserObject & pps = static_cast<const NodalOverlapUserObject &>(y);
}
