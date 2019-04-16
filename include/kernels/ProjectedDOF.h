//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef PROJECTEDDOF_H
#define PROJECTEDDOF_H

#include "NodalKernel.h"
#include "NodalOverlapUserObject.h"
#include "NodalDOFUserObject.h"

class ProjectedDOF;

template <>
InputParameters validParams<ProjectedDOF>();

class ProjectedDOF : public NodalKernel
{
public:
    ProjectedDOF(const InputParameters & parameters);

    //Compute the residual function
    virtual Real computeQpResidual() override;

protected:

    unsigned int n_macro_dof = 12; //!The number of degrees of freedom for each macro node
    unsigned int n_DNS_dof = 3; //!The number of degrees of freedom for each micro node

    //A user object which contains the map between node number and order in the DOF vector
    const NodalOverlapUserObject& _nodal_overlap;

    //A user object which contains the DOF vectors from the DNS
    const NodalDOFUserObject& _DNS_dof;

    //A user object which contains the DOF vectors from the macro-scale
    const NodalDOFUserObject& _macro_dof;
    
    //The degree of freedom the kernel applies to
    unsigned int _dof_num;

    //A flag indicating whether the nodes are associated with the DNS or macro-scale
    bool _is_DNS;

    //A scale factor multipled by the DOF values to ensure the residual has a non-negligible size
    double _scale_factor;

    //Maps from the node number to the order in the projection matrices
    const std::map< dof_id_type, unsigned int >* macro_node_to_col;
    const std::map< dof_id_type, unsigned int >* micro_node_to_row;

    const overlap::EigVec* Dh_micro; //The macro-scale ghost dof vector obtained from DNS dof projection
    //const overlap::EigVec* Dh_macro; //The macro-scale ghost dof vector obtained from macro dof projection //All zero so ignored for now
    const overlap::EigVec* Qh_micro; //The DNS ghost dof vector obtained from DNS dof projection
    const overlap::EigVec* Qh_macro; //The macro-scale ghost dof vector obtained from DNS dof projection
};

#endif /* PROJECTEDDOF_H */
