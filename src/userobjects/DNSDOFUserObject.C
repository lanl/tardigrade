//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DNSDOFUserObject.h"

registerMooseObject("tardigradeApp", DNSDOFUserObject);

template <>
InputParameters
validParams<DNSDOFUserObject>()
{
    InputParameters params = validParams<NodalUserObject>();
    params.addRequiredCoupledVar("u1", "The displacement in the 1 direction");
    params.addRequiredCoupledVar("u2", "The displacement in the 2 direction");
    params.addRequiredCoupledVar("u3", "The displacement in the 3 direction");
    params.addRequiredParam<UserObjectName>("nodal_overlap_userobject", "The userobject associated with the nodal overlap.");
    params.addRequiredParam<UserObjectName>("projector_userobject", "The userobject which defines the projection between lengthscales.");
    params.addParam<UserObjectName>("micromorphic_DOF_userobject", "The userobject which collects the micromorphic degrees of freedom.");
    return params;
}

DNSDOFUserObject::DNSDOFUserObject(const InputParameters & parameters)
    : NodalUserObject(parameters),
    _u1(coupledValue("u1")),
    _u2(coupledValue("u2")),
    _u3(coupledValue("u3")),
    _nodal_overlap(getUserObject<NodalOverlapUserObject>("nodal_overlap_userobject")),
    _projector(getUserObject<ProjectorUserObject>("projector_userobject")),
    _micromorphic_DOF(getUserObject<MicromorphicDOFUserObject>("micromorphic_DOF_userobject"))
{
}

void
DNSDOFUserObject::initialize()
{
    _console << "Initializing DNS DOF UserObject: " << name() << std::endl;

   //Get the macro to row and micro to col maps
   micro_node_to_row = _nodal_overlap.get_micro_node_to_row(); 

   //Get the nodal information
   _nodal_overlap.get_node_info(num_macro_ghost, num_macro_free, num_micro_ghost, num_micro_free);

   //Initialize the free degree of freedom vector
   Q = std::vector< double >(num_micro_dof*num_micro_free, 0);//overlap::EigVec::Zero(num_micro_dof*num_micro_free);

   //Initialize the ghost degree of freedom vectors
   Dh = std::vector< double >(num_macro_dof*num_macro_ghost, 0);//overlap::EigVec::Zero(num_macro_dof*num_macro_ghost);
   Qh = std::vector< double >(num_micro_dof*num_micro_ghost, 0);//overlap::EigVec::Zero(num_micro_dof*num_micro_ghost);
   
}

void
DNSDOFUserObject::execute()
{

    std::map< dof_id_type, unsigned int >::const_iterator it;

    //Check if the current node is in the micro domain
    it = micro_node_to_row->find(_current_node->id());
    if (it != micro_node_to_row->end()){
        // If the node is ``free'' perform projection
        if (it->second<num_micro_free){

            Q[num_micro_dof*it->second + 0] = _u1[0];
            Q[num_micro_dof*it->second + 1] = _u2[0];
            Q[num_micro_dof*it->second + 2] = _u3[0];
        }
        return;
    }
    return;
}

void
DNSDOFUserObject::threadJoin(const UserObject & y)
{
}

void
DNSDOFUserObject::finalize()
{

    std::cout << "Finalizing DNS DOF\n";
    const std::vector< double >* D = _micromorphic_DOF.get_D();
    _projector.project_dof((*D), Q, Dh, Qh);

//    std::cout << "Q:\n";
//    for (unsigned int i=0; i<Q.size()/3; i++){
//        if (fabs(Q[3*i+1])>1e-4){
//            std::cout << "index, value: " << 3*i+1 << ", " << Q[3*i+1] << "\n";
//        }
//    }

//    //Get the solver object
//    BDhQsolver = _projector.get_BDhQsolver();
//
//    //Get the shapefunction matrix
//    const overlap::SpMat* shapefunction = _projector.get_shapefunction();
//
//    //Extract the sub-shapefunction matrices
//    //Get the micro-information
//    unsigned int num_macro_free, num_macro_ghost, num_micro_free, num_micro_ghost;
//    _nodal_overlap.get_node_info(num_macro_ghost, num_macro_free, num_micro_ghost, num_micro_free);
//
//    overlap::SpMat NQhD = shapefunction->block( num_micro_dof*num_micro_free, 0, num_micro_dof*num_micro_ghost,  num_macro_dof*num_macro_free);
//    overlap::SpMat NQhDh = shapefunction->block( num_micro_dof*num_micro_free, num_macro_dof*num_macro_free, num_micro_dof*num_micro_ghost, num_macro_dof*num_macro_ghost);
//
//    //Solve for Dh
//    Dh = BDhQsolver->solve(Q);
//
//    std::cout << "Dh:\n" << Dh << "\n";
//
//    //Solve for Qh
//    const overlap::EigVec* D_macro = _micromorphic_DOF.get_D();
//    Qh = NQhD*(*D_macro) + NQhDh*Dh;
//
//    std::cout << "Qh:\n" << Qh << "\n";
//
    std::cout << "End of Nodal DOF\n\n";
}

const std::vector< double >* DNSDOFUserObject::get_Dh() const{
    /*!
    Return a pointer to the ghost macro dof vector
    */

    return &Dh;
}

const std::vector< double >* DNSDOFUserObject::get_Qh() const{
    /*!
    Return a pointer to the ghost DNS dof vector
    */

    return &Qh;
}
