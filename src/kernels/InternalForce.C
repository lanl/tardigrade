/*!
=====================================================================
|                         InternalForce.cpp                         |
---------------------------------------------------------------------
| The source file for the micromorphic internal force kernel. This  |
| kernel draws upon the balance equations from the git repository   |
| micromorphic_element available at:                                |
|     https://bitbucket.org/NateAM2/                                |
| and requires twelve degrees of freedom to be defined for solid    |
| mechanics. These degrees of freedom are:                          |
|     u_1, u_2, u_3, phi_11, phi_22, phi_33, phi_23, phi_13, phi,12 |
|     phi_32, phi_31, phi_21                                        |
| where u_1 -> u_3 are macro deformations, and phi_11 -> phi_21 are |
| the components of the micro-displacement tensor.                  |
=====================================================================
*/

#include<InternalForce.h>

//We define the valid parameters for this kernel and their default values
registerMooseObject("tardigradeApp", InternalForce);

template<>
InputParameters
validParams<InternalForce>(){
    InputParameters params = validParams<Kernel>();
    params.addRequiredParam<int>("component", "The component of the internal force vector");
    params.addRequiredParam<int>("dof_num",   "The degree of freedom to use for the diagonal jacobian calculation");
    return params;
}

InternalForce::InternalForce(const InputParameters & parameters)
    : // We have to call the constructor for the base class first
        Kernel(parameters),
        _component(getParam<int>("component")),
        _dof_num(getParam<int>("dof_num")),
        _cauchy(getMaterialProperty<std::vector<double>>("cauchy")),
        _DcauchyDgrad_u(getMaterialProperty<std::vector<std::vector<double>>>("DcauchyDgrad_u")),
        _DcauchyDphi(getMaterialProperty<std::vector<std::vector<double>>>("DcauchyDphi")),
        _DcauchyDgrad_phi(getMaterialProperty<std::vector<std::vector<double>>>("DcauchyDgrad_phi"))
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

Real InternalForce::computeQpResidual(){
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
    
    //Copy the test function so that the balance equation function can read it
    double dNdx[3];
    for (int indx=0; indx<3; indx++){dNdx[indx] = _grad_test[_i][_qp](indx);}//p+i);}

    balance_equations::compute_internal_force(_component, dNdx, _cauchy[_qp], fint);
    //std::cout << "fint: " << fint << "\n";
    return fint;
}

Real InternalForce::computeQpJacobian(){
    /*!
    ===========================
    |    computeQpJacobian    |
    ===========================

    Compute the diagonal jacobian term.

    */
    Real dfdUint;

    //Copy the test and interpolation functions so that the balance equation function can read it
    double dNdx[3];
    double detadx[3];
    for (int indx=0; indx<3; indx++){
        dNdx[indx]   = _grad_test[_i][_qp](indx);
        detadx[indx] = _grad_phi[_i][_qp](indx);
    }

    balance_equations::compute_internal_force_jacobian(_component,           _dof_num, 
                                                       _test[_i][_qp],       dNdx, _phi[_j][_qp],          detadx,
                                                       _DcauchyDgrad_u[_qp], _DcauchyDphi[_qp],   _DcauchyDgrad_phi[_qp], dfdUint);

//    std::cout << "N:    " << _test[_i][_qp] << "\n";
//    std::cout << "dNdx: ";
//    for (int indx=0; indx<3; indx++){ std::cout << dNdx[indx] << " ";}
//    std::cout << "\n";
//
//    std::cout << "eta:    " << _phi[_i][_qp] << "\n";
//    std::cout << "detadx: ";
//    for (int indx=0; indx<3; indx++){ std::cout << detadx[indx] << " ";}
//    std::cout << "\n";
//
//    std::cout << "_DcauchyDgrad_u:\n";
//    for (int indx=0; indx<_DcauchyDgrad_u[_qp].size(); indx++){
//        for (int jndx=0; jndx<_DcauchyDgrad_u[_qp][indx].size(); jndx++){
//            std::cout << _DcauchyDgrad_u[_qp][indx][jndx] << " ";
//        }
//        std::cout << "\n";
//    }
//
//    std::cout << "_DcauchyDphi:\n";
//    for (int indx=0; indx<_DcauchyDphi[_qp].size(); indx++){
//        for (int jndx=0; jndx<_DcauchyDphi[_qp][indx].size(); jndx++){
//            std::cout << _DcauchyDphi[_qp][indx][jndx] << " ";
//        }
//        std::cout << "\n";
//    }
//
//    std::cout << "_DcauchyDgrad_phi:\n";
//    for (int indx=0; indx<_DcauchyDgrad_phi[_qp].size(); indx++){
//        for (int jndx=0; jndx<_DcauchyDgrad_phi[_qp][indx].size(); jndx++){
//            std::cout << _DcauchyDgrad_phi[_qp][indx][jndx] << " ";
//        }
//        std::cout << "\n";
//    }
//
//    std::cout << "dfdUint: " << dfdUint << "\n";
//    mooseError("fail");
    return dfdUint;
}
