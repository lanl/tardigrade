/*!
====================================================================
|                 MicromorphicLinearElasticity.cpp                 |
====================================================================
| The source file for a class which computes the cauchy stress and |
| the associated jacobians for a linear elastic micromorphic       |
| material.                                                        |
--------------------------------------------------------------------
| Notes: Relies on libraries from the micromorphic_element         |
|        repository available at bitbucket.org/NateAM2             |
====================================================================
*/


#include "MicromorphicLinearElasticity.h"

registerMooseObject("tartigrade", MicromorphicLinearElasticity);

template<>
InputParameters
validParams<MicromorphicLinearElasticity>(){
    InputParameters params = validParams<Material>();

    // Vectors of material properties
    params.addRequiredParam<std::vector<Real>>(
        "material_parameters", "The vector of material parameters required for the stiffness matrices");

    return params;
}

MicromorphicLinearElasticity::MicromorphicLinearElasticity(const InputParameters & parameters)
    : Material(parameters),
    // Declare that this material is going to provide Eigen matrices containing the cauchy stress and 
    // jacobians that Kernels can use.
    _cauchy(declareProperty<Vector_9>("cauchy")),
    _DcauchyDgrad_u(declareProperty<Matrix_9x9>("DcauchyDgrad_u")),
    _DcauchyDphi(declareProperty<Matrix_9x9>("DcauchyDphi")),
    _DcauchyDgrad_phi(declareProperty<Matrix_9x27>("DcauchyDgrad_phi"))
    _parameters(getParam<std::vector<Real>>("material_parameters"){}

void MicromorphicLinearElasticity::computeQpProperties(){
    /*!
    =============================
    |    computeQpProperties    |
    =============================

    Evaluate the constitutive model at the quadrature point.

    */

    
}
