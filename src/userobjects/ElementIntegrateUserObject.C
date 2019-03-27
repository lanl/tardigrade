//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementIntegrateUserObject.h"

#include "libmesh/bounding_box.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/quadrature.h"

registerMooseObject("tardigradeApp", ElementIntegrateUserObject);

template <>
InputParameters
validParams<ElementIntegrateUserObject>()
{
    InputParameters params = validParams<ElementUserObject>();
    params.addRequiredCoupledVar("variable", "A dummy variable. Use one of the displacements.");
    params.addParam<MaterialPropertyName>("density", 1., "The name of the property providing density.");
    return params;
}

ElementIntegrateUserObject::ElementIntegrateUserObject(const InputParameters & parameters)
    : ElementUserObject(parameters),
      MooseVariableInterface<Real>(this,
                                   false,
                                   "variable",
                                   Moose::VarKindType::VAR_ANY,
                                   Moose::VarFieldType::VAR_FIELD_STANDARD),
      _var(*mooseVariable()),
      _test(_var.phi()),
      _density(getMaterialProperty<Real>("density"))
{
}

void
ElementIntegrateUserObject::initialize()
{
    _console << "Initializing Element Overlap UserObject: " << name() << std::endl;
}

void
ElementIntegrateUserObject::execute()
{

    computeIntegral();

}

void
ElementIntegrateUserObject::computeIntegral()
{
    std::vector< double > nodal_volumes(_test.size(), 0);
    std::vector< double > nodal_weights(_test.size(), 0);
    std::vector< double > nodal_densities(_test.size(), 0);

    //Integrate the functions
    for (_qp = 0; _qp < _qrule->n_points(); _qp++){
        for (_i = 0; _i < _test.size(); _i++){
            nodal_volumes[_i] += _JxW[_qp]*_coord[_qp]*_test[_i][_qp];
            nodal_weights[_i] += _JxW[_qp]*_coord[_qp]*_test[_i][_qp];
            nodal_densities[_i] += _JxW[_qp]*_coord[_qp]*_test[_i][_qp]*_density[_qp];
        }
    }

    //Add the contributions back to the overall values
    for (_i = 0; _i < _test.size(); _i++){
        std::map< dof_id_type, double>::iterator it;
        it = integrated_weights.find(_current_elem->node_id(_i));

        if (it == integrated_weights.end()){
            integrated_volumes.insert(std::pair< dof_id_type, double >(_current_elem->node_id(_i), nodal_volumes[_i]));
            integrated_weights.insert(std::pair< dof_id_type, double >(_current_elem->node_id(_i), nodal_weights[_i]));
            integrated_weighted_densities.insert(std::pair< dof_id_type, double >(_current_elem->node_id(_i), nodal_densities[_i]));
        }
        else{
            integrated_volumes[_current_elem->node_id(_i)] += nodal_volumes[_i];
            integrated_weights[_current_elem->node_id(_i)] += nodal_weights[_i];
            integrated_weighted_densities[_current_elem->node_id(_i)] += nodal_densities[_i];
        }
    }
    return;
}

int
ElementIntegrateUserObject::get_nodal_volume(dof_id_type micro_node_id, double &volume) const{
    /*!
    Return the nodal volume.
    */

    auto volume_iterator = integrated_volumes.find(micro_node_id);
    
    if (volume_iterator == integrated_volumes.end()){
        mooseError("Node not found");
    }
    else{
        volume = volume_iterator->second;
    }
    return 0;
}

double
ElementIntegrateUserObject::get_nodal_volume(dof_id_type micro_node_id) const{
    /*!
    Return the nodal volume.
    */

    double volume;
    get_nodal_volume(micro_node_id, volume);
    return volume;
}

int
ElementIntegrateUserObject::get_nodal_density(dof_id_type micro_node_id, double &density) const{
    /*!

    Compute the nodal density and return it

    */

    auto weight_iterator = integrated_weights.find(micro_node_id);
    auto density_iterator = integrated_weighted_densities.find(micro_node_id);

    if ((weight_iterator == integrated_weights.end()) || (density_iterator == integrated_weighted_densities.end())){
        mooseError("Node not found");
    }
    else{
        density = density_iterator->second/weight_iterator->second;
    }
    return 0;
}

double
ElementIntegrateUserObject::get_nodal_density(dof_id_type micro_node_id) const{
    /*!
    Compute the nodal density and return it
    */

    double density;
    get_nodal_density(micro_node_id, density);
    return density;
}

Real
ElementIntegrateUserObject::computeQpIntegral()
{
    return 1;
}

void
ElementIntegrateUserObject::threadJoin(const UserObject & y)
{

}

void
ElementIntegrateUserObject::finalize()
{
    std::cout << "Finalizing Element Overlap\n";
//    for (unsigned int i=0; i<integrated_weights.size(); i++){
//        std::cout << integrated_weights[i] << ", " << integrated_weighted_densities[i] << ", " << integrated_weighted_densities[i]/integrated_weights[i]  << "\n";
//    }
//    mooseError("derp");
}
