#include<ProjectedDOF.h>

registerMooseObject("tardigradeApp", ProjectedDOF);

//Define the valid parameters for this kernel and their default values

template<>
InputParameters
validParams<ProjectedDOF>(){
    InputParameters params = validParams<NodalKernel>();
    params.addRequiredParam<UserObjectName>("nodal_overlap_userobject","The UserObject that defines the overlapped nodes");
    params.addRequiredParam<UserObjectName>("DNS_dof_userobject","The userobject that computes the ghost DOF (type DNSDOFUserObject)");
    params.addRequiredParam<unsigned int>("dof_num", "The degree of freedom number (0-2) for displacement (3-11) for micro-displacement");
    params.addParam<bool>("is_DNS", false, "A flag indicating whether the nodes are associated with the DNS or macro-scale");
    params.addParam<double>("scale_factor", 1., "A scale factor to help with convergence if required");
    return params;
}

ProjectedDOF::ProjectedDOF(const InputParameters & parameters)
    : // Call the constructor for the base class
        NodalKernel(parameters),
        _nodal_overlap(getUserObject<NodalOverlapUserObject>("nodal_overlap_userobject")),
        _dof_object(getUserObject<DNSDOFUserObject>("DNS_dof_userobject")),
        _dof_num(getParam<unsigned int>("dof_num")),
        _is_DNS(getParam<bool>("is_DNS")),
        _scale_factor(getParam<double>("scale_factor")){

    /*!
    The constructor for the ProjectedDOF kernel. This takes the degrees of freedom as projected using 
    the macro and DNS NodalDOFUserObjects and applies them via a kernel.

    Note that for this kernel no attempt is (currently) made to implement the Jacobian.
    */

}

Real ProjectedDOF::computeQpResidual(){
    /*!
    Compute the residual at the node for the indicated component
    */

    //Flag for if the result is a DNS
    std::map< dof_id_type, unsigned int>::const_iterator it;
    unsigned int index, num_macro_ghost, num_macro_free, num_micro_ghost, num_micro_free;

    //Get the degree of freedom maps
    macro_node_to_col = _nodal_overlap.get_macro_node_to_col();
    micro_node_to_row = _nodal_overlap.get_micro_node_to_row();

    //Get the information about the number of free and ghost nodes
    _nodal_overlap.get_node_info(num_macro_ghost, num_macro_free, num_micro_ghost, num_micro_free);

    //Get the degrees of freedom
    Dh = _dof_object.get_Dh();
//    Dh_macro = _macro_dof.get_Dh(); //will always be zero for now
    Qh = _dof_object.get_Qh();

//    std::cout << "Dh_micro.size(): " << Dh_micro->size() << "\n";
//    std::cout << "Qh_micro.size(): " << Qh_micro->size() << "\n";
//    std::cout << "Qh_macro.size(): " << Qh_macro->size() << "\n";


    if (_is_DNS){
        
        it = micro_node_to_row->find(_current_node->id());

        //Check if the node is a ghost DNS node in the overlap domain
        if ((it != micro_node_to_row->end()) && (it->second >= num_micro_free)){

            index = it->second - num_micro_free;
            index *= n_DNS_dof;
            index += _dof_num;
//            if (_dof_num == 2){mooseError("A micro dof of 2 has been observed\n");}
            if (index >= Qh->size()){
                return 0;
                mooseError("Error: index " + std::to_string(index) + " is larger than allowed\n" +
                           "Check your inputs for dof_num (currently " + std::to_string(_dof_num) + ")\n" +
                           "to ensure that they are between 0-2 for the DNS terms\n" + 
                           "  Kernel name: " + name() + "\n" +
                           "  node id: " + std::to_string(_current_node->id()) + "\n" + 
                           "  node order: " + std::to_string(it->second) + "\n" + 
                           "  # of DOF associated with DNS node: " + std::to_string(n_DNS_dof) + "\n" + 
                           "  Qh.size(): " + std::to_string(Qh->size()) + "\n");
            }
//            std::cout << "dof_num, Qh_micro, Qh_macro, Qh_micro + Qh_macro: " << _dof_num << ", " << (*Qh_micro)[index] << ", " << (*Qh_macro)[index] << ", " <<_scale_factor*((*Qh_micro)[index] + (*Qh_macro)[index]) << "\n";
            return _scale_factor*(_u[_qp]-(*Qh)[index]);
        }
        return 0;
    }
    else{

        it = macro_node_to_col->find(_current_node->id());

        //Check if the node is a ghost micromorphic node in the overlap domain
        if ((it != macro_node_to_col->end()) && (it->second >= num_macro_free)){

            index = it->second - num_macro_free;
            index *= n_macro_dof;
            index += _dof_num;
//            if (_dof_num == 2){mooseError("A macro dof of 2 has been observed\n");}
            if (index >= Dh->size()){
                return 0;
                mooseError("Error: index " + std::to_string(index) + " is larger than allowed\n" +
                           "Check your inputs for dof_num (currently " + std::to_string(_dof_num) + ")\n" +
                           "to ensure that they are between 0-11 for the micromorphic terms\n" + 
                           "  Kernel name: " + name() + "\n" +
                           "  node id: " + std::to_string(_current_node->id()) + "\n" +
                           "  node order: " + std::to_string(it->second) + "\n" +  
                           "  # of DOF associated with micromorphic node: " + std::to_string(n_macro_dof) + "\n" + 
                           "  Dh.size(): " + std::to_string(Dh->size()) + "\n");
            }
//            std::cout << "node id, dof_num, _u.size(), _u[_qp], Dh, R: " << _current_node->id() << ", " << _dof_num << ", " << _u.size() << ", " << _u[_qp] << ", " << (*Dh)[index] << ", " << _scale_factor*(_u[_qp] - (*Dh)[index]) << "\n";
            return _scale_factor*(_u[_qp]-(*Dh)[index]);
        }
        return 0;
    }

    return 0;
}

Real ProjectedDOF::computeQpJacobian(){

    //Flag for if the result is a DNS
    std::map< dof_id_type, unsigned int>::const_iterator it;
    unsigned int num_macro_ghost, num_macro_free, num_micro_ghost, num_micro_free;

    //Get the degree of freedom maps
    macro_node_to_col = _nodal_overlap.get_macro_node_to_col();
    micro_node_to_row = _nodal_overlap.get_micro_node_to_row();

    //Get the information about the number of free and ghost nodes
    _nodal_overlap.get_node_info(num_macro_ghost, num_macro_free, num_micro_ghost, num_micro_free);

    if (_is_DNS){
        it = micro_node_to_row->find(_current_node->id());

        if ((it != micro_node_to_row->end()) && (it->second >= num_micro_free)){
            return _scale_factor;
        }
    }
    else{
        it = macro_node_to_col->find(_current_node->id());

        if ((it != macro_node_to_col->end()) && (it->second >= num_macro_free)){
            return _scale_factor;
        }
    }

    return 0;

}
