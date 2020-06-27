#pragma once

#include "AuxKernel.h"
#include "ElementIntegrateUserObject.h"

class NodalVolumeAverages;

template <>
InputParameters validParams<NodalVolumeAverages>();

class NodalVolumeAverages : public AuxKernel
{

    public:
        NodalVolumeAverages( const InputParameters & parameters );

    protected:
        virtual Real computeValue( ) override;

        const ElementIntegrateUserObject & _elementIntegrate;

        const bool & _compute_nodal_volume;

};
