/*=============================================================================
 *                      OverlapCouplingGeneralUserObject                      *
 *============================================================================*
 * Communicate with the overlap coupling libraries outside of MOOSE to        *
 * perform the overlap coupling. The general schme is the following:          *
 *   - Prior to first increment                                               *
 *     - OverlapCouplingNodeUserObject: Initialize overlap variables to empty *
 *                                      ( both macro and micro )              *
 *   - Timestep Begins                                                        *
 *     - OverlapCouplingNodeUserObject: Initialize overlap variables to       *
 *                                      empty ( both macro and micro )        *
 *     - Nonlinear iteration                                                  *
 *       - Linear iteration( s )                                              *
 *         - All balance equation kernels update                              *
 *         - CoupledPenaltyDirichletBC: Apply nodal BC values from the        *
 *                                      OverlapCouplingNodeUserObject         *
 *                                      ( both macro and micro )              *
 *                                                                            *
 *       - OverlapCouplingGeneralUserObject: Run the external overlap coupling*
 *                                           method ( macro only ) and update *
 *                                           the degree of freedom maps ( both*
 *                                           macro and micro )                *
 *============================================================================*/

#ifndef OVERLAPCOUPLINGNODALUSEROBJECT
#define OVERLAPCOUPLINGNODALUSEROBJECT
#include "GeneralUserObject.h"

class OverlapCouplingGeneralUserObject;

template<>
InputParameters validParams<GeneralUserObject>( );

class OverlapCouplingGeneralUserObject : public GeneralUserObject {

    public:
        OverlapCouplingGeneralUserObject( const InputParameters & parameters );

        virtual void initialize( ) override;

        virtual void execute( ) override;

        virtual void threadJoin( const UserObject & y ) override;

        virtual void finalize( ) override;

        //!Return the value
        bool getNodalValue( const dof_id_type nodeId, const dof_id_type dofIndex, Real &nodalDOFValue );

    protected:

        const std::string &_overlap_coupling_filename;

        unsigned int _dim;

        const bool &_isMacro;

        std::vector< dof_id_type > _nodeIds;

        std::vector< double > _dofValues;

};

#endif
