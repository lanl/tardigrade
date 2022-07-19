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

InputParameters
OverlapCoupling::validParams()
{
    InputParameters params = GeneralUserObject::validParams( );
    params.addRequiredParam< bool > ( "is_macroscale", "Flag for whether the coupling object is in the macroscale or not" );
    params.addParam< std::string >( "overlap_configuration_filename", "_none_", "The overlap configuration YAML filename" );
    return params;

}

OverlapCoupling::OverlapCoupling( const InputParameters &parameters )
    : GeneralUserObject( parameters ),
        _isMacro( getParam< bool > ( "is_macroscale" ) ),
        _overlapCouplingFilename( getParam< std::string >( "overlap_configuration_filename" ) )
{

    if ( _isMacro ){

        if ( _overlapCouplingFilename.compare( "_none_" ) == 0 ){

            mooseError( "The 'overlap_configuration_filename' MUST be defined for the macroscale" );

        }

    }
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

    _console << "map size: " << _microGlobalLocalNodeMap.size( ) << "\n";
    _console << "dof size: " << _updatedMicroDisplacementDOF.size( ) << "\n";

    return;

}

void OverlapCoupling::execute( ){
    /*!
     * Execute the general user object
     */

    std::cout << "execute\n";

    if ( !_isMacro ){
        _console << name( ) << " is not macro\n";
        return;
    }
    _console << name( ) << " is macro\n";

    mooseInfo( "Performing overlap coupling" );

    overlapCoupling::errorOut error = overlapCoupling::runOverlapCoupling( _overlapCouplingFilename,
                                                                           _microGlobalLocalNodeMap, _updatedMicroDisplacementDOF,
                                                                           _microAugmentedLagrangianForce,
                                                                           _macroGlobalLocalNodeMap, _updatedMacroDisplacementDOF,
                                                                           _macroAugmentedLagrangianForce );

    if ( error ){

        error->print( );

        delete error;

        mooseError( "failure in performing the overlap coupling" );

    }

    _console << "_microGlobalLocalNodeMap.size( ):     " << _microGlobalLocalNodeMap.size( ) << "\n";
    _console << "_updatedMicroDisplacementDOF.size( ): " << _updatedMicroDisplacementDOF.size( ) << "\n";

    std::cerr << "exiting execute\n";

//    mooseError( "The overlap coupling configuration file should be updated" );

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

void OverlapCoupling::setAttribute( const std::unordered_map< unsigned int, unsigned int > &attribute, const std::string &attributeName ){
    /*!
     * Set the attribute of the class
     */

    if ( attributeName.compare( "microGlobalLocalNodeMap" ) == 0 ){

        _microGlobalLocalNodeMap = attribute;

    }
    else{

        mooseError( "Attribute name " + attributeName + " not recognized." );

    }
}

void OverlapCoupling::setAttribute( const std::vector< double > &attribute, const std::string &attributeName ){
    /*!
     * Set the attribute of the class
     */

    if ( attributeName.compare( "updatedMicroDisplacementDOF" ) == 0 ){

        _updatedMicroDisplacementDOF = attribute;

    }
    else if ( attributeName.compare( "microAugmentedLagrangianForce" ) == 0 ){
        
        _microAugmentedLagrangianForce = attribute;

    }
    else{

        mooseError( "Attribute name " + attributeName + " not recognized." );

    }
}

const std::unordered_map< unsigned int, unsigned int > *OverlapCoupling::getMicroGlobalLocalNodeMap( ) const{
    /*!
     * Return a reference to the micro dof map
     */

    return &_microGlobalLocalNodeMap;
}

const std::vector< double > *OverlapCoupling::getMicroDisplacementDOF( ) const{
    /*!
     * Return a reference to the micro displacement
     */

    return &_updatedMicroDisplacementDOF;
}

const std::unordered_map< unsigned int, unsigned int > *OverlapCoupling::getMacroGlobalLocalNodeMap( ) const{
    /*!
     * Return a reference to the macro dof map
     */

    return &_macroGlobalLocalNodeMap;
}

const std::vector< double > *OverlapCoupling::getMacroDisplacementDOF( ) const{
    /*!
     * Return a reference to the macro displacement
     */

    return &_updatedMacroDisplacementDOF;
}

const std::vector< double > *OverlapCoupling::getMacroAugmentedLagrangianForce( ) const{
    /*!
     * Return a reference to the macro augmented lagrangian force
     */

    return &_macroAugmentedLagrangianForce;
}

const std::vector< double > *OverlapCoupling::getMicroAugmentedLagrangianForce( ) const{
    /*!
     * Return a reference to the micro augmented lagrangian force
     */

    return &_microAugmentedLagrangianForce;
}
