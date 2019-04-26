/*!
=====================================================================
|                     InternalCoupleOverlap.h                       |
---------------------------------------------------------------------
| The header file for the micromorphic internal couple kernel when  |
| used in an overlap coupling environment. This kernel draws upon   |
| the balance equations from the git repository                     |
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


#ifndef INTERNALCOUPLEOVERLAP_H
#define INTERNALCOUPLEOVERLAP_H

#include "InternalCouple.h"
#include <balance_equations.h> //The Micromorphic balance equations
#include "NodalOverlapUserObject.h"

//Forward declarations
class InternalCoupleOverlap;

template <>
InputParameters validParams<InternalCoupleOverlap>();

class InternalCoupleOverlap : public InternalCouple{
    /*!
    -------------------------------------------------
    |    InternalCoupleOverlap kernel definition    |
    -------------------------------------------------

    We define the kernel for the internal couple to be

    -(psi (PK2_{ij} - SIGMA_{ij}) - psi_{,k} M_{kji})

    where psi is the test function, PK2 is the second piola-kirchoff stress,
    SIGMA is the symmetric stress, and M is the higher order couple stress.
    Note that this is a vector definition so we will need to define multiple 
    kernels each obtaining a specific component of the residual.
    */
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        //The constructor definition
        InternalCoupleOverlap(const InputParameters & parameters);

    protected:
        //The residual at a quadrature point.
        //Note: As stated above, a kernel is a scalar quantity so we have to define
        //      as many of these as dimensions in our problem.
        virtual Real computeQpResidual() override;
        
        //The diagonal members of the jacobian at a quadrature point.
        //Note: The parameter chosen as the "diagonal" member is totally arbitrary
        //      We just need to have each DOF associated with some residual term.
        virtual Real computeQpJacobian() override;

        //The off-diagonal members of the jacobian at the quadrature point.
        virtual Real computeQpOffDiagJacobian(unsigned int jvar) override; 

        //A user object which contains the map between node number and order in the DOF vector
        const NodalOverlapUserObject& _nodal_overlap;

        //A function to determine if the current node is a ghost node or not
        bool node_is_ghost();
};

#endif
