/*!
=====================================================================
|                     InternalCoupleOverlap.h                       |
---------------------------------------------------------------------
| The header file for the micromorphic internal couple kernel used  |
| in an overlap-coupling environment. This kernel draws upon the    |
| balance equations from the git repository  micromorphic_element   |
| available at:                                                     |
|     https://bitbucket.org/NateAM2/                                |
| and requires twelve degrees of freedom to be defined for solid    |
| mechanics. These degrees of freedom are:                          |
|     u_1, u_2, u_3, phi_11, phi_22, phi_33, phi_23, phi_13, phi,12 |
|     phi_32, phi_31, phi_21                                        |
| where u_1 -> u_3 are macro deformations, and phi_11 -> phi_21 are |
| the components of the micro-displacement tensor.                  |
=====================================================================
*/

#include<InternalCoupleOverlap.h>

//We define the valid parameters for this kernel and their default values
registerMooseObject("tardigradeApp", InternalCoupleOverlap);

template<>
InputParameters
validParams<InternalCoupleOverlap>(){
//    std::cout << "In InputParameters InternalCoupleOverlap\n";
    InputParameters params = validParams<InternalCouple>();
    params.addRequiredParam<UserObjectName>("nodal_overlap_userobject","The UserObject that defines the overlapped nodes");
    return params;
}

InternalCoupleOverlap::InternalCoupleOverlap(const InputParameters & parameters)
    : // We have to call the constructor for the base class first
        InternalCouple(parameters),
        _nodal_overlap(getUserObject<NodalOverlapUserObject>("nodal_overlap_userobject"))
    {
    /*!
    =====================
    |    Constructor    |
    =====================

    The constructor for the InternalForce class.
    Note that this constructor is just variable 
    assignments.

    */
}

Real InternalCoupleOverlap::computeQpResidual(){
    /*!
    ===========================
    |    computeQpResidual    |
    ===========================

    Compute the residual at the quadrature point for 
    the indicated component.

    cint_ij = psi (PK2_ij - SIGMA_ij) -psi_k M_kji

    where i = _component_i
          j = _component_j

    */
    Real cint;
    if (node_is_ghost()){
        return 0;
    }
    cint = InternalCouple::computeQpResidual();
    return cint;
}

Real InternalCoupleOverlap::computeQpJacobian(){
    /*!
    ===========================
    |    computeQpJacobian    |
    ===========================

    Compute the diagonal jacobian term.

    */
    Real dcdUint;
    if (node_is_ghost()){
        return 0;
    }
    dcdUint = InternalCouple::computeQpJacobian();
    return dcdUint;
}

Real InternalCoupleOverlap::computeQpOffDiagJacobian(unsigned int jvar){
    /*!
    ==================================
    |    computeQpOffDiagJacobian    |
    ==================================

    Compute the off-diagonal terms of the jacobian
    */


    Real dcdUint;
    if (node_is_ghost()){
        return 0;
    }
    dcdUint = InternalCouple::computeQpOffDiagJacobian(jvar);
    return dcdUint;
}

bool InternalCoupleOverlap::node_is_ghost(){
    /*
    Check if the current node is a ghost node
    */

    const Node*_node = _current_elem->node_ptr(_i); //Get a pointer to the node associated with test function i

    //Check if the node is in the overlap region
    const std::map< dof_id_type, unsigned int>* macro_node_to_col = _nodal_overlap.get_macro_node_to_col();
    std::map< dof_id_type, unsigned int>::const_iterator it = macro_node_to_col->find(_node->id());

    if (it == macro_node_to_col->end()){
        return false;
    }

    //Check if the node is ghost
    unsigned int num_macro_ghost, num_macro_free, num_micro_ghost, num_micro_free;
    _nodal_overlap.get_node_info(num_macro_ghost, num_macro_free, num_micro_ghost, num_micro_free);

    if (it->second < num_macro_free){
        return false;
    }

    //The node is ghost
    return true;
}

