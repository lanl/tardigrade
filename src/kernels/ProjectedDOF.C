#include<ProjectedDOF.h>

registerMooseObject("tardigradeApp", ProjectedDOF);

//Define the valid parameters for this kernel and their default values

template<>
InputParameters
validParams<ProjectedDOF>(){
    InputParameters params = validParams<NodalKernel>();
    params.addRequiredParam<UserObjectName>("nodal_overlap_userobject","The UserObject that defines the overlapped nodes");
    params.addRequiredParam<UserObjectName>("DNS_dof_userobject","The nodal degree of freedom object associated with the DNS");
    params.addRequiredParam<UserObjectName>("macro_dof_userobject","The nodal degree of freedom object associated with the macro-scale");
    params.addRequiredParam<unsigned int>("dof_num", "The degree of freedom number (0-2) for displacement (3-11) for micro-displacement");
    params.addParam<bool>("is_DNS", false, "A flag indicating whether the nodes are associated with the DNS or macro-scale");
    params.addParam<double>("scale_factor", 1., "A scale factor to help with convergence if required");
    return params;
}

ProjectedDOF::ProjectedDOF(const InputParameters & parameters)
    : // Call the constructor for the base class
        NodalKernel(parameters),
        _nodal_overlap(getUserObject<NodalOverlapUserObject>("nodal_overlap_userobject")),
        _DNS_dof(getUserObject<NodalDOFUserObject>("DNS_dof_userobject")),
        _macro_dof(getUserObject<NodalDOFUserObject>("macro_dof_userobject")),
        _dof_num(getParam<unsigned int>("dof_num")),
        _is_DNS(getParam<bool>("is_DNS")),
        _scale_factor(getParam<double>("scale_factor")){

    /*!
    The constructor for the ProjectedDOF kernel. This takes the degrees of freedom as projected using 
    the macro and DNS NodalDOFUserObjects and applies them via a kernel.

    Note that for this kernel no attempt is (currently) made to implement the Jacobian.
    */

    //Get the degree of freedom maps
    macro_node_to_col = _nodal_overlap.get_macro_node_to_col();
    micro_node_to_row = _nodal_overlap.get_micro_node_to_row();

    //Get the degrees of freedom
    Dh_micro = _DNS_dof.get_Dh();
//    Dh_macro = _macro_dof.get_Dh(); //will always be zero for now
    Qh_micro = _DNS_dof.get_Qh();
    Qh_macro = _macro_dof.get_Qh();
}

Real ProjectedDOF::computeQpResidual(){
    /*!
    Compute the residual at the node for the indicated component
    */

    //Flag for if the result is a DNS
    unsigned int index = _current_node->id();
    if (_is_DNS){
        index *= n_DNS_dof;
        index += _dof_num;
        return _scale_factor*((*Qh_micro)[index] + (*Qh_macro)[index]);
    }
    else{
        index *= n_macro_dof;
        index += _dof_num;
        return _scale_factor*(*Dh_micro)[index];
    }
}
