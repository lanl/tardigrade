//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef NODALOVERLAPUSEROBJECT_H
#define NODALOVERLAPUSEROBJECT_H

// MOOSE includes
#include "NodalUserObject.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/bounding_box.h"
#include "libmesh/fe_interface.h"

// Forward Declarations
class NodalOverlapUserObject;
class ContainedNode;

template <>
InputParameters validParams<NodalOverlapUserObject>();

class NodalOverlapUserObject : public NodalUserObject{
public:
  NodalOverlapUserObject(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override;

protected:

    //!Parameters
    SubdomainName _macroscale_domain; //! The name of the macro-scale subdomain

//    MooseMesh _projection_mesh; //! The mesh to be used in the projection
    SubdomainID _macro_id; //! The id of the macro-scale subdomain
    BoundingBox _macro_bounding_box; //! The bounding box of the macro-scale subdomain
    std::map<dof_id_type, std::vector< dof_id_type > > macro_element_to_micro_nodes; //! A map from the macro-scale element id to the contained micro-scale node ids
    std::vector< dof_id_type > micro_elements; //! A vector of elements which are connected to the contained nodes.
};

#endif

class ContainedNode{
    /*!
    A class which defines the required values for a node which is found to lie within an 
    element.
    */

    public:

    libMesh::Node* node;
    libMesh::Point local_coordinates;

    ContainedNode(){}

    ContainedNode(libMesh::Node* _node, libMesh::Point _local_coordinates)
    {
        node = _node;
        local_coordinates = _local_coordinates;
    }

};
