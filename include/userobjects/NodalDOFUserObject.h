//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef NODALDOFUSEROBJECT_H
#define NODALDOFUSEROBJECT_H

// MOOSE includes
#include "NodalUserObject.h"
#include "NodalOverlapUserObject.h"
#include "ProjectorUserObject.h"

// Forward Declarations
class NodalDOFUserObject;
class ContainedNode;

template <>
InputParameters validParams<NodalDOFUserObject>();

class NodalDOFUserObject : public NodalUserObject{
    public:
        NodalDOFUserObject(const InputParameters & parameters);

        virtual void initialize() override;
        virtual void execute() override;
        virtual void threadJoin(const UserObject & y) override;
        virtual void finalize() override;

        const overlap::EigVec* get_Dh() const;
        const overlap::EigVec* get_Qh() const;

    protected:

        //!The number of micro degrees of freedom
        unsigned int num_micro_dof = 3;

        //!The number of macro degrees of freedom
        unsigned int num_macro_dof = 12;

        //!The variables to project
        const VariableValue& _u1;
        const VariableValue& _u2;
        const VariableValue& _u3;
        const VariableValue& _phi11;
        const VariableValue& _phi22;
        const VariableValue& _phi33;
        const VariableValue& _phi23;
        const VariableValue& _phi13;
        const VariableValue& _phi12;
        const VariableValue& _phi32;
        const VariableValue& _phi31;
        const VariableValue& _phi21;

        //!Flag indicating if the domain is the DNS
        const bool& _is_DNS;

        //!The nodal overlap userobject
        const NodalOverlapUserObject & _nodal_overlap;

        //!The projector userobject
        const ProjectorUserObject & _projector;

        //Projection information
        unsigned int num_macro_ghost;
        unsigned int num_macro_free;
        unsigned int num_micro_ghost;
        unsigned int num_micro_free;

        //Projection Matrices
        const overlap::SpMat* BDhD;
        const overlap::SpMat* BDhQ;
        const overlap::SpMat* BQhD;
        const overlap::SpMat* BQhQ;

        //Maps from the node number to the order in the projection matrices
        const std::map< dof_id_type, unsigned int >* macro_node_to_col;
        const std::map< dof_id_type, unsigned int >* micro_node_to_row;

        //The computed contributions to the degree of freedom vectors
        overlap::EigVec Dh;
        overlap::EigVec Qh;

};

#endif
