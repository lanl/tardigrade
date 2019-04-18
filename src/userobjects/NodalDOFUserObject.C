//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NodalDOFUserObject.h"

registerMooseObject("tardigradeApp", NodalDOFUserObject);

template <>
InputParameters
validParams<NodalDOFUserObject>()
{
    InputParameters params = validParams<NodalUserObject>();
    params.addRequiredCoupledVar("u1", "The displacement in the 1 direction");
    params.addRequiredCoupledVar("u2", "The displacement in the 2 direction");
    params.addRequiredCoupledVar("u3", "The displacement in the 3 direction");
    params.addCoupledVar("phi11", 0, "The micro-deformation in the 11 direction");
    params.addCoupledVar("phi22", 0, "The micro-deformation in the 22 direction");
    params.addCoupledVar("phi33", 0, "The micro-deformation in the 33 direction");
    params.addCoupledVar("phi23", 0, "The micro-deformation in the 23 direction");
    params.addCoupledVar("phi13", 0, "The micro-deformation in the 13 direction");
    params.addCoupledVar("phi12", 0, "The micro-deformation in the 12 direction");
    params.addCoupledVar("phi32", 0, "The micro-deformation in the 32 direction");
    params.addCoupledVar("phi31", 0, "The micro-deformation in the 31 direction");
    params.addCoupledVar("phi21", 0, "The micro-deformation in the 21 direction");
    params.addRequiredParam<bool>("is_DNS", "Indicates if the node is in the DNS or not");
    params.addRequiredParam<UserObjectName>("nodal_overlap_userobject", "The userobject associated with the nodal overlap.");
    params.addRequiredParam<UserObjectName>("projector_userobject", "The userobject which defines the projection between lengthscales.");
    return params;
}

NodalDOFUserObject::NodalDOFUserObject(const InputParameters & parameters)
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
    _is_DNS(getParam<bool>("is_DNS")),
    _nodal_overlap(getUserObject<NodalOverlapUserObject>("nodal_overlap_userobject")),
    _projector(getUserObject<ProjectorUserObject>("projector_userobject"))
{
}

void
NodalDOFUserObject::initialize()
{
    _console << "Initializing Nodal DOF UserObject: " << name() << std::endl;

   //Get the macro to row and micro to col maps
   macro_node_to_col = _nodal_overlap.get_macro_node_to_col();
   micro_node_to_row = _nodal_overlap.get_micro_node_to_row(); 

   //Get the nodal information
   _nodal_overlap.get_node_info(num_macro_ghost, num_macro_free, num_micro_ghost, num_micro_free);

   //Get the projectors
   BDhD = _projector.get_BDhD();
   BDhQ = _projector.get_BDhQ();
   BQhD = _projector.get_BQhD();
   BQhQ = _projector.get_BQhQ();

   //Initialize the degree of freedom vectors
   Dh = overlap::EigVec::Zero(num_macro_dof*num_macro_ghost);
   Qh = overlap::EigVec::Zero(num_micro_dof*num_micro_ghost);
   
}

void
NodalDOFUserObject::execute()
{

    std::map< dof_id_type, unsigned int >::const_iterator it;

    //Check if the current node is in the micro domain
    if (_is_DNS){

        //Check if the current node is in the micro domain
        it = micro_node_to_row->find(_current_node->id());
        if (it != micro_node_to_row->end()){
            // If the node is ``free'' perform projection
            if (it->second<num_micro_free){
                overlap::SpEigVec Qp = overlap::SpEigVec(num_micro_dof*num_micro_free);
                Qp.insert(num_micro_dof*it->second + 0) = _u1[0];
                Qp.insert(num_micro_dof*it->second + 1) = _u2[0];
                Qp.insert(num_micro_dof*it->second + 2) = _u3[0];

//                std::cout << BQhD->row(num_micro_dof*it->second + 1) << "\n";
//                mooseError("NARF!");

                //Add the contributions to the DOF vectors
                overlap::SpEigVec res = (*BDhQ)*Qp;
                Dh += res;
                std::cout << "_u2: " << _u2[0] << "\n";
                res = (*BQhQ)*Qp;
                Qh += res;
//                Dh += (*BDhQ)*Qp;
//                Qh += (*BQhQ)*Qp;
            }
            return;
        }
    }
    //Check if the current node is in the macro domain
    else {
        it = macro_node_to_col->find(_current_node->id());
        if (it != macro_node_to_col->end()){

            // If the node is ``free'' perform projection
            if (it->second<num_macro_free){
                overlap::SpEigVec Dp = overlap::SpEigVec(num_macro_dof*num_macro_free);
                Dp.insert(num_macro_dof*it->second +  0) = _u1[0];
                Dp.insert(num_macro_dof*it->second +  1) = _u2[0];
                Dp.insert(num_macro_dof*it->second +  2) = _u3[0];
                Dp.insert(num_macro_dof*it->second +  3) = _phi11[0];
                Dp.insert(num_macro_dof*it->second +  4) = _phi22[0];
                Dp.insert(num_macro_dof*it->second +  5) = _phi33[0];
                Dp.insert(num_macro_dof*it->second +  6) = _phi23[0];
                Dp.insert(num_macro_dof*it->second +  7) = _phi13[0];
                Dp.insert(num_macro_dof*it->second +  8) = _phi12[0];
                Dp.insert(num_macro_dof*it->second +  9) = _phi32[0];
                Dp.insert(num_macro_dof*it->second + 10) = _phi31[0];
                Dp.insert(num_macro_dof*it->second + 11) = _phi21[0];


                //Add the contributions to the DOF vectors (Ignore BDhatD for now)
                overlap::SpEigVec res = (*BQhD)*Dp;
                Qh += res;
//                Qh = (*BQhD)*Dp;
//                if (fabs(_u3[0])>1e-5){
//                    std::cout << "Dp:\n" << Dp << "\n";
//                    std::cout << "Qh:\n" << Qh << "\n";
//                    mooseError("Narf!");
//                }
            }
    
            return;
        }
    }
/*    if (fabs(_u2[0]) > 0.005){
        std::cout << "_u1: "; for (unsigned int i=0; i<_u1.size(); i++){std::cout << _u1[i] << " ";} std::cout << "\n";
        std::cout << "_u2: "; for (unsigned int i=0; i<_u2.size(); i++){std::cout << _u2[i] << " ";} std::cout << "\n";
        std::cout << "_u3: "; for (unsigned int i=0; i<_u3.size(); i++){std::cout << _u3[i] << " ";} std::cout << "\n";
        std::cout << "Qh:\n" << Qh << "\n";
        mooseError("Narf!");
    }
*/    return;
}

void
NodalDOFUserObject::threadJoin(const UserObject & y)
{
}

void
NodalDOFUserObject::finalize()
{

    std::cout << "Finalizing Nodal DOF\n";

//    std::cout << "Qh:\n" << Qh << "\n";
//    std::cout << "Dh:\n" << Dh << "\n";

    std::cout << "End of Nodal DOF\n\n";
}

const overlap::EigVec* NodalDOFUserObject::get_Dh() const{
    /*!
    Return a pointer to the ghost macro dof vector
    */

    return &Dh;
}

const overlap::EigVec* NodalDOFUserObject::get_Qh() const{
    /*!
    Return a pointer to the ghost DNS dof vector
    */

    return &Qh;
}
