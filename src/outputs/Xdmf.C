/*!
 * Output the results of the simulation in XDMF format
 */

#include "Xdmf.h"

registerMooseObject( "tardigradeApp", Xdmf );

template<>
InputParameters
validParams<Xdmf>(){
    InputParameters params = validParams<AdvancedOutput>();
    params.addClassDescription( "Output of system in XDMF format" );
//    params += AdvancedOutput::enableOutputTypes( "" );

    return params;
}

Xdmf::Xdmf( const InputParameters &parameters )
    : AdvancedOutput( parameters ){

    std::cout << "in constructor\n";

    return;
}

void Xdmf::initialSetup(){

    std::cout << "in initial setup\n";

//    shared_ptr< XdmfDomain >  _root = XdmfDomain::New();
//    std::cout << "forming writer\n";
//    shared_ptr< XdmfWriter > writer = XdmfWriter::New( filename() );
//    std::cout << "accepting writer\n";
//    _root->accept( writer );

}
