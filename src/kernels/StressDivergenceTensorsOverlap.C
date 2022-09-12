/*!
=====================================================================
|                     StressDivergenceTensorsOverlap.cpp            |
---------------------------------------------------------------------
| The source file for the micromorphic coupling capable Stress      |
| Divergence Kernel. This class inherits from the                   |
| StressDivergenceTensors kernel.                                   |
=====================================================================
*/

#include<StressDivergenceTensorsOverlap.h>

//We define the valid parameters for this kernel and their default values
registerMooseObject("tardigradeApp", StressDivergenceTensorsOverlap);

template<>
InputParameters
validParams<StressDivergenceTensorsOverlap>(){
    InputParameters params = validParams<StressDivergenceTensors>();
    params.addRequiredParam<UserObjectName>("nodal_overlap_userobject","The UserObject that defines the overlapped nodes");
    return params;
}

StressDivergenceTensorsOverlap::StressDivergenceTensorsOverlap(const InputParameters & parameters)
    : // We have to call the constructor for the base class first
        StressDivergenceTensors(parameters),
        _nodal_overlap(getUserObject<NodalOverlapUserObject>("nodal_overlap_userobject"))
    {
    /*!
    =====================
    |    Constructor    |
    =====================

    The constructor for the StressDivergenceTensorsOverlap class.
    Note that this constructor is just variable 
    assignments.

    */

}

Real StressDivergenceTensorsOverlap::computeQpResidual(){
    /*!
    ===========================
    |    computeQpResidual    |
    ===========================

    Compute the residual at the quadrature point for 
    the indicated component.

    fint_i = -psi_j sigma_ji

    where i = _component

    */

    Real fint;

    //Check if the current node is ghost
    if (node_is_ghost()){
        //The node is ghost. There is no stress residual.
        return 0;
    }

    fint = StressDivergenceTensors::computeQpResidual();
    return fint;
}

Real StressDivergenceTensorsOverlap::computeQpJacobian(){
    /*!
    ===========================
    |    computeQpJacobian    |
    ===========================

    Compute the diagonal jacobian term.

    */

    Real dfdUint;

    //Check if the current node is ghost
    if (node_is_ghost()){
        return 0;
    }
    dfdUint = StressDivergenceTensors::computeQpJacobian();
    return dfdUint;
}

Real StressDivergenceTensorsOverlap::computeQpOffDiagJacobian(unsigned int jvar){
    /*!
    ==================================
    |    computeQpOffDiagJacobian    |
    ==================================

    Compute the off-diagonal terms of the jacobian
    */

    Real dfdUint;
    //Check if the current node is ghost
    if (node_is_ghost()){
        return 0;
    }
    dfdUint = StressDivergenceTensors::computeQpOffDiagJacobian(jvar);
    return dfdUint;
}

bool StressDivergenceTensorsOverlap::node_is_ghost(){
    /*
    Check if the current node is a ghost node
    */

    const Node*_node = _current_elem->node_ptr(_i); //Get a pointer to the node associated with test function i

    //Check if the node is in the overlap region
    const std::map< dof_id_type, unsigned int>* micro_node_to_row = _nodal_overlap.get_micro_node_to_row();
    std::map< dof_id_type, unsigned int>::const_iterator it = micro_node_to_row->find(_node->id());

    if (it == micro_node_to_row->end()){
        return false;
    }

    //Check if the node is ghost
    unsigned int num_macro_ghost, num_macro_free, num_micro_ghost, num_micro_free;
    _nodal_overlap.get_node_info(num_macro_ghost, num_macro_free, num_micro_ghost, num_micro_free);

    if (it->second < num_micro_free){
        return false;
    }

    //The node is ghost
    return true;
}
