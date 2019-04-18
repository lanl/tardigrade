//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MicromorphicDOFUserObject.h"

registerMooseObject("tardigradeApp", MicromorphicDOFUserObject);

template <>
InputParameters
validParams<MicromorphicDOFUserObject>()
{
    InputParameters params = validParams<NodalUserObject>();
    params.addRequiredCoupledVar("u1", "The displacement in the 1 direction");
    params.addRequiredCoupledVar("u2", "The displacement in the 2 direction");
    params.addRequiredCoupledVar("u3", "The displacement in the 3 direction");
    params.addRequiredCoupledVar("phi11", "The micro-deformation in the 11 direction");
    params.addRequiredCoupledVar("phi22", "The micro-deformation in the 22 direction");
    params.addRequiredCoupledVar("phi33", "The micro-deformation in the 33 direction");
    params.addRequiredCoupledVar("phi23", "The micro-deformation in the 23 direction");
    params.addRequiredCoupledVar("phi13", "The micro-deformation in the 13 direction");
    params.addRequiredCoupledVar("phi12", "The micro-deformation in the 12 direction");
    params.addRequiredCoupledVar("phi32", "The micro-deformation in the 32 direction");
    params.addRequiredCoupledVar("phi31", "The micro-deformation in the 31 direction");
    params.addRequiredCoupledVar("phi21", "The micro-deformation in the 21 direction");
    params.addRequiredParam<UserObjectName>("nodal_overlap_userobject", "The userobject associated with the nodal overlap.");
    params.addRequiredParam<UserObjectName>("projector_userobject", "The userobject which defines the projection between lengthscales.");
    return params;
}

MicromorphicDOFUserObject::MicromorphicDOFUserObject(const InputParameters & parameters)
    : NodalUserObject(parameters),
    _u1(coupledValue("u1")),
    _u2(coupledValue("u2")),
    _u3(coupledValue("u3")),
    _phi11(coupledValue("phi11")),
    _phi22(coupledValue("phi22")),
    _phi33(coupledValue("phi33")),
    _phi23(coupledValue("phi23")),
    _phi13(coupledValue("phi13")),
    _phi12(coupledValue("phi12")),
    _phi32(coupledValue("phi32")),
    _phi31(coupledValue("phi31")),
    _phi21(coupledValue("phi21")),
    _nodal_overlap(getUserObject<NodalOverlapUserObject>("nodal_overlap_userobject")),
    _projector(getUserObject<ProjectorUserObject>("projector_userobject")),
{
}

void
MicromorphicDOFUserObject::initialize()
{
    _console << "Initializing Micromorphic DOF UserObject: " << name() << std::endl;

   //Get the macro to row and micro to col maps
   macro_node_to_col = _nodal_overlap.get_macro_node_to_col();

   //Get the nodal information
   _nodal_overlap.get_node_info(num_macro_ghost, num_macro_free, num_micro_ghost, num_micro_free);

   //Initialize the free degree of freedom vector
   D = overlap::EigVec::Zero(num_macro_dof*num_macro_free);

   //Initialize the ghost degree of freedom vectors
}

void
MicromorphicDOFUserObject::execute()
{

    std::map< dof_id_type, unsigned int >::const_iterator it;

    it = macro_node_to_col->find(_current_node->id());
    if (it != macro_node_to_col->end()){

        // If the node is ``free'' perform projection
        if (it->second<num_macro_free){

            D[num_macro_dof*it->second +  0] = _u1[0];
            D[num_macro_dof*it->second +  1] = _u2[0];
            D[num_macro_dof*it->second +  2] = _u3[0];
            D[num_macro_dof*it->second +  3] = _phi11[0];
            D[num_macro_dof*it->second +  4] = _phi22[0];
            D[num_macro_dof*it->second +  5] = _phi33[0];
            D[num_macro_dof*it->second +  6] = _phi23[0];
            D[num_macro_dof*it->second +  7] = _phi13[0];
            D[num_macro_dof*it->second +  8] = _phi12[0];
            D[num_macro_dof*it->second +  9] = _phi32[0];
            D[num_macro_dof*it->second + 10] = _phi31[0];
            D[num_macro_dof*it->second + 11] = _phi21[0];
        }
    
        return;
    }
    return;
}

void
MicromorphicDOFUserObject::threadJoin(const UserObject & y)
{
}

void
MicromorphicDOFUserObject::finalize()
{

    std::cout << "Finalizing Micromorphic DOF\n";

    std::cout << "End of Micromorphic DOF\n\n";
}

const overlap::EigVec* MicromorphicDOFUserObject::get_D() const{
    /*!
    Return a pointer to the free macro dof vector
    */

    return &D;
}
