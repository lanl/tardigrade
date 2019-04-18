//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef DNSDOFUSEROBJECT_H
#define DNSDOFUSEROBJECT_H

// MOOSE includes
#include "NodalUserObject.h"
#include "NodalOverlapUserObject.h"
#include "ProjectorUserObject.h"
#include "MicromorphicDOFUserObject.h"

// Forward Declarations
class DNSDOFUserObject;
class ContainedNode;

template <>
InputParameters validParams<DNSDOFUserObject>();

class DNSDOFUserObject : public NodalUserObject{
    public:
        DNSDOFUserObject(const InputParameters & parameters);

        virtual void initialize() override;
        virtual void execute() override;
        virtual void threadJoin(const UserObject & y) override;
        virtual void finalize() override;

        const overlap::EigVec* get_Dh() const;
        const overlap::EigVec* get_Qh() const;

    protected:

        //!The number of macro degrees of freedom
        unsigned int num_macro_dof = 12;

        //!The number of micro degrees of freedom
        unsigned int num_micro_dof = 3;

        //!The variables to project
        const VariableValue& _u1;
        const VariableValue& _u2;
        const VariableValue& _u3;

        //!The nodal overlap userobject
        const NodalOverlapUserObject & _nodal_overlap;

        //!The projector userobject
        const ProjectorUserObject & _projector;

        //!The micromorphic DOF userobject
        const MicromorphicDOFUserObject & _micromorphic_DOF;

        //Projection information
        unsigned int num_macro_ghost;
        unsigned int num_macro_free;
        unsigned int num_micro_ghost;
        unsigned int num_micro_free;

        //Projection Matrices and solvers
        const overlap::QRsolver* BDhQsolver;

        //Maps from the node number to the order in the projection matrices
        const std::map< dof_id_type, unsigned int >* micro_node_to_row;

        //The degree of freedom vectors
        overlap::EigVec Q; 

        //The computed contributions to the degree of freedom vectors
        overlap::EigVec Dh;
        overlap::EigVec Qh;

};

#endif
