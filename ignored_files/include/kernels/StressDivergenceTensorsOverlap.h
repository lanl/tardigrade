/*!
=====================================================================
|                 StressDivergenceTensorsOverlap.h                  |
---------------------------------------------------------------------
| The header file for the StressDivergenceTensors kernel which has  |
| been extended to work in a micromoprhic overlap coupling problem. |
=====================================================================
*/


#ifndef STRESSDIVERGENCETENSORSOVERLAP_H
#define STRESSDIVERGENCETENSORSOVERLAP_H

#include "StressDivergenceTensors.h"
#include "NodalOverlapUserObject.h"

//Forward declarations
class StressDivergenceTensorsOverlap;

template <>
InputParameters validParams<StressDivergenceTensorsOverlap>();

class StressDivergenceTensorsOverlap : public StressDivergenceTensors{
    /*!
    ----------------------------------------------------------
    |    StressDivergenceTensorsOverlap kernel definition    |
    ----------------------------------------------------------

    */
    public:
        //The constructor definition
        StressDivergenceTensorsOverlap(const InputParameters & parameters);

    protected:
        //The residual at a quadrature point.
        //Note: As stated above, a kernel is a scalar quantity so we have to define
        //      as many of these as dimensions in our problem.
        virtual Real computeQpResidual() override;
        
        //The diagonal members of the residual at a quadrature point.
        //Note: The parameter chosen as the "diagonal" member is totally arbitrary
        //      We just need to have each DOF associated with some residual term.
        virtual Real computeQpJacobian() override;

        //The off diagonal members of the residual at a quadrature point.
        virtual Real computeQpOffDiagJacobian(unsigned int jvar) override; 

        //A user object which contains the map between node number and order in the DOF vector
        const NodalOverlapUserObject& _nodal_overlap;

        //A function to determine if a node is a ghost node or not
        bool node_is_ghost();
};

#endif
