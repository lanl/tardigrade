//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef PROJECTORUSEROBJECT_H
#define PROJECTORUSEROBJECT_H

// MOOSE includes
#include "NodalOverlapUserObject.h"
#include "ElementIntegrateUserObject.h"
#include "occonfiguration.h"
#include "overlap_coupling.h"

// Forward Declarations
class ProjectorUserObject;

template <>
InputParameters validParams<ProjectorUserObject>();

class ProjectorUserObject : public ElementUserObject{
    public:
        ProjectorUserObject(const InputParameters & parameters);

        virtual void initialize() override;
        virtual void execute() override;
        virtual void threadJoin(const UserObject & y) override;
        virtual void finalize() override;

        const overlap::SpMat* get_shapefunction() const;
        const overlap::QRsolver* get_BDhQsolver() const;
        const overlap::QRsolver* get_BDhQ_transpose_solver() const;

    protected:
        //Settings
        unsigned int n_macro_dof = 12; //!The number of degrees of freedom for each macro node
        unsigned int n_micro_dof = 3; //!The number of degrees of freedom for each micro node
        bool solve_for_projectors = true; //!Whether to solve for the projection matrices or not

        //Required user objects
        const NodalOverlapUserObject & _nodal_overlap;
        const ElementIntegrateUserObject & _element_integrate;

        //Variables required for the computation of the projector
        std::map< dof_id_type, std::vector< double > > volumes; //! The volume of the DNS in each gauss domain (reference config)
        std::map< dof_id_type, std::vector< double > > densities; //! The density of the DNS in each gauss domain (reference config)
        std::map< dof_id_type, overlap::vecOfvec > cgs; //! The centers of gravity of the DNS in each gauss domain (reference config)

        //Projections from the macro space to the micro space
        const std::map< dof_id_type, unsigned int >* macro_node_to_col;
        const std::map< dof_id_type, unsigned int >* micro_node_to_row;

        //Map from nodes located on element boundaries to the number of elements they should be incorporated in the shape functions
        const std::map< dof_id_type, unsigned int >* micro_node_elcount;

        //Define the dns weights objects
        std::map < dof_id_type, std::vector< overlap::integrateMap > > dns_weights;
        
        //Define the vectors which will be used to store the entries to the shapefunction matrix
        std::vector< overlap::T > tripletList;

        //Define the shapefunction matrix
        overlap::SpMat shapefunction;

        //Define the projectors
        overlap::QRsolver BDhQsolver;
        overlap::QRsolver BDhQ_transpose_solver;

        //!Utility Methods
        void collect_local_nodes(overlap::vecOfvec &local_nodes, std::vector< dof_id_type > &macro_node_ids);
        void collect_quadrature_points(overlap::vecOfvec &local_gpts);
        void collect_local_micronode_positions(overlap::vecOfvec &local_micronode_positions);
        void map_integrateMap_values( unsigned int gpt);
        void compute_dns_weights(const overlap::vecOfvec &local_nodes, const overlap::vecOfvec &local_gpts,
                                 const std::vector< dof_id_type >* micro_nodes, const overlap::vecOfvec &local_micronode_positions);
        void get_cg_local_coordinates(unsigned int gpt, Point &local_cg);
        void get_cg_phi(unsigned int gpt, overlap::vecOfvec &phi);
        void compute_shapefunction_matrix_contributions(const std::vector< dof_id_type > &macro_node_ids);

        template< class myType >
        void print_vector(const std::vector< myType > &);
        template< class myType >
        void print_matrix(const std::vector< std::vector< myType > > &);
        void solve_for_projector(const overlap::SpMat &A, const overlap::SpMat &B, overlap::SpMat &X);

};

#endif
