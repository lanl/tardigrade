#pragma once

#include "AuxKernel.h"
#include "ElementIntegrateUserObject.h"

class NodalVolumeAverages;

class NodalVolumeAverages : public AuxKernel
{

    public:
        NodalVolumeAverages( const InputParameters & parameters );
        
        static InputParameters validParams();

    protected:
        virtual Real computeValue( ) override;

        const ElementIntegrateUserObject & _elementIntegrate;

        const bool & _compute_nodal_volume;

};
