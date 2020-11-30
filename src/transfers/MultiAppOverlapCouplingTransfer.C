
#include "MultiAppOverlapCouplingTransfer.h"
#include "UserObject.h"
#include "MultiApp.h"

registerMooseObject( "tardigradeApp", MultiAppOverlapCouplingTransfer );

defineLegacyParams( MultiAppTransfer );

template< >
InputParameters validParams< MultiAppOverlapCouplingTransfer >( ){

    InputParameters params = MultiAppTransfer::validParams( );

    params.addRequiredParam< UserObjectName >( "source_user_object", "The OverlapCoupling object you want to transfer values from" );
    params.addRequiredParam< UserObjectName >( "target_user_object", "The OverlapCoupling object you want to transfer values to" );

    return params;

}

MultiAppOverlapCouplingTransfer::MultiAppOverlapCouplingTransfer( const InputParameters & parameters )
    : MultiAppTransfer( parameters ),
      _source_user_object_name( getParam< UserObjectName >( "source_user_object" ) ),
      _target_user_object_name( getParam< UserObjectName >( "target_user_object" ) ){

        _fe_problem.mesh( ).errorIfDistributedMesh( "MultiAppOverlapCouplingTransfer" );

        return;

    }

void MultiAppOverlapCouplingTransfer::execute( ){

    _console << "Beginning MultiAppOverlapCouplingTransfer " << name( ) << std::endl;

    switch ( _current_direction ){

            case TO_MULTIAPP: {
                for ( unsigned int i = 0; i < _multi_app->numGlobalApps( ); i++ ){

                    Moose::ScopedCommSwapper swapper( _multi_app->comm( ) );

                    //Get the OverlapCoupling user objects
                    const OverlapCoupling &main_user_object =
                        _multi_app->problemBase( ).getUserObject< OverlapCoupling >( _source_user_object_name );

                    _multi_app->appProblemBase( i ).getUserObject< OverlapCoupling >( _target_user_object_name ).setAttribute( *main_user_object.getMicroGlobalLocalNodeMap( ), "microGlobalLocalNodeMap" );
                    _multi_app->appProblemBase( i ).getUserObject< OverlapCoupling >( _target_user_object_name ).setAttribute( *main_user_object.getMicroDisplacementDOF( ), "updatedMicroDisplacementDOF" );
                    _multi_app->appProblemBase( i ).getUserObject< OverlapCoupling >( _target_user_object_name ).setAttribute( *main_user_object.getMicroAugmentedLagrangianForce( ), "microAugmentedLagrangianForce" );

                }

            }
//
//            case FROM_MULTIAPP:
//                    mooseError( "FROM_MULTIAPP not supported" );


    }

    _console << "MultiAppOverlapCouplingTransfer " << name ( ) << " completed." << std::endl;

    return;

}
