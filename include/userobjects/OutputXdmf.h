/*!
 * Output the results of the simulation in XDMF format
 */

#ifndef OUTPUTXDMF_HEADER
#define OUTPUTXDMF_HEADER

//MOOSE includes
#include "GeneralUserObject.h"
#include "MooseMesh.h"
#include "AdvancedOutput.h"

//XDMF includes
#include "XdmfUnstructuredGrid.hpp"

//Forward declarations
class OutputXdmf;

/**
 * Class for output data to the XDMF format
 */

class OutputXdmf : public GeneralUserObject
{
    public:
        static InputParameters validParams( );

        /*!
         * Class constructor
         */

        OutputXdmf( const InputParameters & parameters );

        /*!
         * Initialization method
         */

        virtual void initialize( ) override;

        /*!
         * Override the execute method
         */

        virtual void execute( ) override;

        /*!
         * Override the finalize method
         */

        virtual void finalize( ) override;

        /*!
         * Output the data
         */

        void output( const ExecFlagType &type );

        /*!
         * Detect of the mesh has changed
         */

        virtual void meshChanged( );

        /*!
         * Get the names of the nodal variables to be output
         */

        const std::set< std::string > &getNodalVariableOutput( );

        /*!
         * Class destructor
         */

        ~OutputXdmf(){};

    protected:
        //Protected  attributes
//        shared_ptr< XdmfDomain > _domain;
//        std::shared_ptr<XdmfHDF5Writer> heavyWriter;
        const std::string &_filename;

        Real _time;

        MooseMesh &_mesh;

        EquationSystems &_es;

        void writeMeshToFile( const bool local = false );

        void outputNodalVariables( const std::set< std::string > * system_names = nullptr );

    private:
        //Private attributes
        bool _xdmf_mesh_changed;
        unsigned int _num_temporal_collections = 0;
        unsigned int _current_reference_grid = 0;
        OutputDataWarehouse _execute_data;

};

#endif
