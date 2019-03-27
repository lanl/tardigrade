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

    //!Get the list of micro-nodes contained in the current macro element
    const std::vector< dof_id_type>* get_relevant_micro_nodes(dof_id_type macro_elem_id) const;

    //!Get the local position of the micro-node in the macro element
    const std::vector< Point >* get_local_node_positions(dof_id_type macro_element_id) const;

protected:

    //!Parameters
    SubdomainName _macroscale_domain; //! The name of the macro-scale subdomain

//    MooseMesh _projection_mesh; //! The mesh to be used in the projection
    //! The id of the macro-scale subdomain
    SubdomainID _macro_id;

     //! The bounding box of the macro-scale subdomain
    BoundingBox _macro_bounding_box;

    //! A map from the macro-scale element id to the contained micro-scale node ids
    std::map<dof_id_type, std::vector< dof_id_type > > macro_element_to_micro_nodes;

    //! A map from the macro-scale element id to the local positions of the contained micro-scale nodes
    std::map<dof_id_type, std::vector< Point > > macro_element_to_micro_positions;

     //! A map from the micro-elements which are connected to the contained nodes to the macro-elements they are associated with.
    std::map< dof_id_type, std::vector< dof_id_type > > micro_elements;

    //!Methods

     //!Add an element-node pair to the macro_element_to_micro_nodes map
    void updateMacroMicroMap(dof_id_type macro_elem_id, dof_id_type node_id, Point local_position);

    //!Add micro-elements which haven't been identified to the micro_elements vector
    void updateMicroElements(dof_id_type macro_elem_id, dof_id_type elem_id);

};

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

#endif
