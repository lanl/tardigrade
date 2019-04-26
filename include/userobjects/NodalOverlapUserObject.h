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

        //!Get the map from the micro nodes to the column index
        const std::map< dof_id_type, unsigned int >* get_micro_node_to_row() const;

        //!Get the map from the macro nodes to the row index
        const std::map< dof_id_type, unsigned int >* get_macro_node_to_col() const;

        //!Return quantities of interest about the nodes contained in the overlap domain
        void get_node_info(unsigned int &_num_macro_ghost, unsigned int &_num_macro_free,
                           unsigned int &_num_micro_ghost, unsigned int &_num_micro_free) const;    

    protected:

        //!Settings
        //!These settings may (eventually) be opened up for the user to modify. For now, they are held fixed in code.
        bool restrict_nodes_to_one_element = false; //! If true a micro-node can only appear in one macro-element. This case only occurs when nodes are exactly on the boundary between two elements.
        int ghost_depth = 1; //! The depth at which the macro-scale element (and all its nodes) will be identified as ghost. Note that ghost will supersede free for macro-scale nodes and free will supersede ghost for micro-scale nodes. This is done so that no free micro-scale nodes will be located within the support of a free macro-scale node.
        double _tolr = libMesh::TOLERANCE; //! The relative tolerance for detecting if points are located inside of an element

        //!Parameters
        SubdomainName _macroscale_domain; //! The name of the macro-scale subdomain

//    MooseMesh _projection_mesh; //! The mesh to be used in the projection
        //! The id of the macro-scale subdomain
        SubdomainID _macro_id;

        //! The bounding box of the macro-scale subdomain
        BoundingBox _macro_bounding_box;

        //! A map from the macro-scale element id to the contained micro-scale node ids
        std::map<dof_id_type, std::vector< dof_id_type > > macro_element_to_micro_nodes;

        //! A map from the macro-scale element to its neighbors
        std::map< dof_id_type, std::vector< dof_id_type > > macro_element_neighbors;

        //! A map from the macro-scale element to the id's of its nodes
        std::map< dof_id_type, std::vector< dof_id_type > > macro_element_nodes;

        //! A map from the macro-scale element id to the local positions of the contained micro-scale nodes
        std::map<dof_id_type, std::vector< Point > > macro_element_to_micro_positions;

        //! A map from the micro-elements which are connected to the contained nodes to the macro-elements they are associated with.
        std::map< dof_id_type, std::vector< dof_id_type > > micro_elements;

        //! A map from the micro-scale node to its column index in the shape-function matrix
        std::map< dof_id_type, unsigned int > micro_node_to_row;

        //! A map from the macro-scale node to its column index in the shape-function matrix
        std::map< dof_id_type, unsigned int > macro_node_to_col;

        //! The dimensions required of the shape-function matrix
        unsigned int num_macro_ghost = 0;
        unsigned int num_micro_ghost = 0;
        unsigned int num_macro_free = 0;
        unsigned int num_micro_free = 0;

        //!Methods

        //!Add an element-node pair to the macro_element_to_micro_nodes map
        void updateMacroMicroMap(dof_id_type macro_elem_id, dof_id_type node_id, Point local_position);

        //!Add the neighbors of the macro-element to the macro_element_neighbors map
        void updateMacroNeighborMap( const Elem* macro_element);

        //!Add micro-elements which haven't been identified to the micro_elements vector
        void updateMicroElements(dof_id_type macro_elem_id, dof_id_type elem_id);

        //!Add the number of nodes in the macro element to the count map
        void updateMacroNodeIds(const Elem* macro_element);

        //!Add the id number of the micro-node to the micro_node_to_row map
        void updateMicroNodeRowMap(const dof_id_type micro_node_id);

        //!Add the id numbers of the macro-element's nodes to the macro_node_to_col map
        void updateMacroNodeColMap(const Elem* macro_element);

        //!Check element depth
        void checkElementDepth(const dof_id_type &macro_element, std::map< dof_id_type, int > &element_depths);

        //!Print a vector (debugging)
        template< class myType >
        void print_vector(std::vector< myType > &v);

        //!Print a matrix (debugging)
        template< class myType >
        void print_matrix(std::vector< std::vector< myType > > &m);

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
