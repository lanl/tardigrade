/*!
 * Output the results of the simulation in XDMF format
 */

#ifndef XDMF_HEADER
#define XDMF_HEADER

//MOOSE includes
#include "AdvancedOutput.h"

//XDMF includes
#include "XdmfUnstructuredGrid.hpp"

//Forward declarations
class Xdmf;

template <>
InputParameters validParams<Xdmf>();

/**
 * Class for output data to the XDMF format
 */

class Xdmf : public AdvancedOutput
{
    public:
        static InputParameters validParams();

        /*!
         * Class constructor
         */

        Xdmf( const InputParameters & parameters );

        /*!
         * Initialization method
         */

        void initialSetup();

        /*!
         * Output the data
         */

        void output( const ExecFlagType &type );

        /*!
         * Class destructor
         */

        ~Xdmf(){};

    protected:
        //Protected  attributes
//        shared_ptr< XdmfDomain > _domain;
//        std::shared_ptr<XdmfHDF5Writer> heavyWriter;
        std::string _filename;

        void writeMeshToFile( const bool local = false );

        void outputNodalVariables( const std::set< std::string > * system_names = nullptr );

};

#endif
