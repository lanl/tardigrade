//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DNSFVectors.h"

#include "libmesh/bounding_box.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/quadrature.h"

registerMooseObject("tardigradeApp", DNSFVectors);

template <>
InputParameters
validParams<DNSFVectors>()
{
    InputParameters params = validParams<ElementUserObject>();
    params.addRequiredCoupledVar("variable", "A dummy variable. Use one of the displacements.");
params.addRequiredParam<UserObjectName>("nodal_overlap_userobject", "The NodalOverlapUserObject that detects which micro-scale nodes overlap with the macro-scale");
    params.addParam<std::string>("base_name", "Material property base name");
    return params;
}

DNSFVectors::DNSFVectors(const InputParameters & parameters)
    : ElementUserObject(parameters),
      MooseVariableInterface<Real>(this,
                                   false,
                                   "variable",
                                   Moose::VarKindType::VAR_ANY,
                                   Moose::VarFieldType::VAR_FIELD_STANDARD),
      _JxW(_assembly.JxW()),
      _coord(_assembly.coordTransformation()),
      _var(*mooseVariable()),
      _test(_var.phi()),
      _grad_test(_var.gradPhi()),
      _nodal_overlap(getUserObject<NodalOverlapUserObject>("nodal_overlap_userobject")),
      _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
      _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress"))
{
}

void
DNSFVectors::initialize()
{
    _console << "Initializing DNSFVectors UserObject: " << name() << std::endl;

    //Get the micro to row maps
    micro_node_to_row = _nodal_overlap.get_micro_node_to_row();

    //Get the DNS information
    _nodal_overlap.get_node_info(num_macro_ghost, num_macro_free, num_micro_ghost, num_micro_free);

    //Size the FintQh vector and initialize it to zero
    FintQh = std::vector< double > (num_micro_ghost*num_micro_dof, 0);
}

void
DNSFVectors::execute()
{

    computeFVectors();

}

void
DNSFVectors::computeFVectors()
{
    /*!
    Compute the F vectors
    */
    std::vector< double > Element_Fint(num_micro_dof*_test.size(), 0);

    dof_id_type node_id;

    //Integrate the functions
    for (_i = 0; _i < _test.size(); _i++){
        //Check if the current node id is ghost
        node_id = _current_elem->node_ptr(_i)->id();
        auto it = micro_node_to_row->find(node_id);

        if ((it != micro_node_to_row->end()) && (it->second >= num_micro_ghost)){

            for (_qp = 0; _qp < _qrule->n_points(); _qp++){
                for (_j = 0; _j < _mesh.dimension(); _j++){
                    Element_Fint[_i+_j] += _JxW[_qp]*_coord[_qp]*computeFintQp();
                }
            }
        }

        //Update residual vectors
        for (_j = 0; _j < _mesh.dimension(); _j++){
            //Update the internal force vector for the ghost DNS DOF
            FintQh[num_micro_dof*it->second + _j] += Element_Fint[_i+_j];
        }

    }

    return;
}

Real
DNSFVectors::computeFintQp()
{
    /*!
    Compute the internal force at the given quadrature point
    */
    return _stress[_qp].row(_j)*_grad_test[_i][_qp];
}

void
DNSFVectors::threadJoin(const UserObject & y)
{

}

void
DNSFVectors::finalize()
{
    std::cout << "Finalizing Element Overlap\n";
//    for (unsigned int i=0; i<integrated_weights.size(); i++){
//        std::cout << integrated_weights[i] << ", " << integrated_weighted_densities[i] << ", " << integrated_weighted_densities[i]/integrated_weights[i]  << "\n";
//    }
//    mooseError("derp");
    std::cout << "End of Element Overlap\n\n";
}
