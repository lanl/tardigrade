/*=============================================================================
 *                               OverlapCoupling                              *
 *============================================================================*
 * Communicate with the overlap coupling libraries outside of MOOSE to        *
 * perform the overlap coupling.                                              *
 *============================================================================*/

#include "OverlapCoupling.h"

#include "overlapCoupling.h"
#include<iostream>
#include<fstream>
#include<chrono>

registerMooseObject( "tardigradeApp", OverlapCoupling );

template<>
InputParameters
validParams< OverlapCoupling >( ){
    InputParameters params = validParams< GeneralUserObject >( );
    params.addRequiredParam< std::string >( "overlap_configuration_filename", "The overlap configuration YAML filename" );
    return params;

}

OverlapCoupling::OverlapCoupling( const InputParameters &parameters )
    : GeneralUserObject( parameters ),
        _overlapCouplingFilename( getParam< std::string >( "overlap_configuration_filename" ) )
{
}

struct cerr_redirect {
    cerr_redirect( std::streambuf * new_buffer )
        : old( std::cerr.rdbuf( new_buffer ) )
    { }

    ~cerr_redirect( ) {
        std::cerr.rdbuf( old );
    }

private:
    std::streambuf * old;
};

void OverlapCoupling::initialize( ){
    /*!
     * Initialize the overlap coupling general user object
     */

    std::cout << "initialize\n";
    std::cout << "execute_enum: " << _current_execute_flag  << "\n";
    _dim = _fe_problem.mesh( ).dimension( );

    return;

}

void OverlapCoupling::execute( ){
    /*!
     * Execute the general user object
     */

    std::cout << "execute\n";

    std::stringbuf buffer;
    cerr_redirect rd( &buffer );

    overlapCoupling::overlapCoupling oc( _overlapCouplingFilename );

    if ( oc.getConstructorError( ) ){

            overlapCoupling::errorOut result
                = new overlapCoupling::errorNode( "OverlapCoupling::execute", "Error in construction of overlapCoupling object" );

        result->addNext( oc.getConstructorError( ) );

        result->print( );

        mooseError( buffer.str( ) );

    }

    overlapCoupling::errorOut error = oc.initializeCoupling( );

    if ( error ){

        overlapCoupling::errorOut result
            = new overlapCoupling::errorNode( "OverlapCoupling::execute", "Error in the initialization of the overlapCoupling object" );

        result->addNext( error );

        result->print( );

        mooseError( buffer.str( ) );

    }

    error = oc.processLastIncrements( );

    if ( error ){

        overlapCoupling::errorOut result
            = new overlapCoupling::errorNode( "OverlapCoupling::execute", "Error in the projection of the data" );

        result->addNext( error );

        result->print( );

        mooseError( buffer.str( ) );

    }

    mooseError( "The overlap coupling configuration file should be updated" );

    return;

}

void OverlapCoupling::threadJoin( const UserObject &y ){
    /*!
     * Join the threads
     */

    std::cout << "thread join\n";

    return;

}

void OverlapCoupling::finalize( ){
    /*!
     * Finalize the user object
     */

    std::cout << "finalizing\n";

    return;

}
