/*!
 * Output the results of the simulation in XDMF format
 */

#include "Xdmf.h"

#include "XdmfDomain.hpp"
#include "XdmfInformation.hpp"
#include "XdmfReader.hpp"
#include "XdmfWriter.hpp"
#include "XdmfHDF5Writer.hpp"
#include "XdmfGeometry.hpp"
#include "XdmfGridCollection.hpp"
#include "XdmfGridCollectionType.hpp"

#include "FEProblem.h"
#include "FileMesh.h"
#include "MooseApp.h"

registerMooseObject( "tardigradeApp", Xdmf );

template<>
InputParameters
validParams<Xdmf>(){
    InputParameters params = validParams<AdvancedOutput>();
    params.addClassDescription( "Output of system in XDMF format" );

    return params;
}

Xdmf::Xdmf( const InputParameters &parameters )
    : AdvancedOutput( parameters ){

    return;
}

void Xdmf::initialSetup(){

    //TODO: Allow for parallel output
    if( this->processor_id() != 0 ){
        return;
    }

    //Call the base class setup method
    AdvancedOutput::initialSetup();

    //Save the filename
    _filename = filename();

    //Open the output file
    shared_ptr< XdmfDomain > _domain = XdmfDomain::New();
    shared_ptr< XdmfInformation > domaininfo = XdmfInformation::New( "Domain", "Primary data structure from a MOOSE simulation" );
    _domain->insert( domaininfo );

    //Create the first temporal collection
    shared_ptr< XdmfGridCollection > _gridHolder = XdmfGridCollection::New( );
    _gridHolder->setType( XdmfGridCollectionType::Temporal( ) );
    shared_ptr< XdmfInformation > _holderInfo = XdmfInformation::New( "Collection_" + std::to_string( _num_temporal_collections ), "The collection of temporal grids" );
    _gridHolder->insert( _holderInfo );
    _domain->insert( _gridHolder );
    _num_temporal_collections++;

    shared_ptr< XdmfHDF5Writer > heavyWriter = XdmfHDF5Writer::New( _filename + ".h5", true );
    heavyWriter->setReleaseData( true );
    shared_ptr< XdmfWriter > writer = XdmfWriter::New( _filename + ".xdmf", heavyWriter );
    _domain->accept( writer );

}

void Xdmf::output( const ExecFlagType &type ){
    /*!
     * Output the simulation data
     *
     * :param const ExecFlagType &type: The type of execution.
     */

    //TODO: Allow for parallel output
    if( this->processor_id() != 0 ){
        return;
    }

    //Write the nodes to the file
    writeMeshToFile( );

    outputNodalVariables( );

    _xdmf_mesh_changed = false;

    return;
}

void Xdmf::meshChanged( ){
    /*!
     * Detect if the mesh has been changed.
     */

    //Run the parent object's meshChanged function
    AdvancedOutput::meshChanged();

    //Alert the object that the mesh has changed
    _xdmf_mesh_changed = true;
}

void Xdmf::writeMeshToFile( const bool local ){
    /*!
     * Write the mesh data to the file
     *
     * :param const bool local: Whether to use the local nodes or not
     */

    /*!==============================
    | Read the existing output file |
    ===============================*/

    //Open the output file
    shared_ptr< XdmfReader > reader = XdmfReader::New();

    //Read in the XDMF domain    
    shared_ptr< XdmfDomain > domain = shared_dynamic_cast< XdmfDomain >( reader->read( _filename + ".xdmf" ) );

    //Construct the output grid
    shared_ptr< XdmfUnstructuredGrid > grid = XdmfUnstructuredGrid::New();

    //Write the timestamp information
    shared_ptr< XdmfTime > untime = XdmfTime::New( _time );
    shared_ptr< XdmfInformation > untimeinfo = XdmfInformation::New( "Time", "This is the current value of the timestep" );
    untime->insert( untimeinfo );
    grid->setTime( untime );

    if ( ( domain->getGridCollection( 0 )->getNumberUnstructuredGrids() > 0 ) && ( !_xdmf_mesh_changed ) ){

        //Point the new grid to the current reference grid's geometry, topology, and ids
        
        shared_ptr< XdmfUnstructuredGrid > reference_grid
            = domain->getGridCollection( 0 )->getUnstructuredGrid( _current_reference_grid );

        //Reference the geometry and topology to the previous grid
        grid->setGeometry( reference_grid->getGeometry( ) );
        grid->setTopology( reference_grid->getTopology( ) );

        //Get the Node and Element ID references
        grid->insert( reference_grid->getAttribute( "NODEID" ) );
        grid->insert( reference_grid->getAttribute( "ELEMID" ) );
//        grid->insert( reference_grid->getAttribute( "SUBDOMAINID" ) );

        //Get the nodesets
        for ( unsigned int i = 0; i < reference_grid->getNumberSets(); i++ ){
            grid->insert( reference_grid->getSet( i ) );
        }

        //Write the grid to the file
        shared_ptr< XdmfGridCollection > collection = domain->getGridCollection( _num_temporal_collections - 1 );

        collection->insert( grid );

        //Construct the file writers
        shared_ptr< XdmfHDF5Writer > heavyWriter = XdmfHDF5Writer::New( _filename + ".h5" );
        heavyWriter->setReleaseData( true );
        shared_ptr< XdmfWriter > writer = XdmfWriter::New( _filename + ".xdmf", heavyWriter );
        writer->setLightDataLimit( 1 );

        //Write out the data
        domain->accept( writer );

        return;
    }

    /*!===============
    | Save the nodes |
    ================*/

    unsigned int _dim = _mesh_ptr->dimension(); //Get the dimension of the mesh
    unsigned int _nNodes = _mesh_ptr->nNodes();  //Determine the number of nodes

    std::vector< unsigned int > _nodeIds, _elementIds, _subdomainIds;
    _nodeIds.reserve( _nNodes ); //Reserve the node id vector
    std::vector< Real > _coordinates;
    _coordinates.reserve( 3 * _nNodes ); //Reserve the coordinate vector

    //Iterate through the nodes
    for ( auto & node : as_range( local ? _mesh_ptr->localNodesBegin() : _mesh_ptr->getMesh().nodes_begin(),
                                  local ? _mesh_ptr->localNodesEnd() : _mesh_ptr->getMesh().nodes_end() ) ){

        //Save the nodal id
        _nodeIds.push_back( node->id() );

        //Save the coordinates of the node
        _coordinates.push_back( (*node)( 0 ) );
        if ( _dim > 1 ){
            _coordinates.push_back( (*node)( 1 ) );
        }
        else{
            _coordinates.push_back( 0. );
        }

        if ( _dim > 2 ){
            _coordinates.push_back( (*node)( 2 ) );
        }
        else{
            _coordinates.push_back( 0. );
        }
    }

    //Save the nodal id information
    shared_ptr<XdmfAttribute> nodeIds = XdmfAttribute::New();
    nodeIds->setType( XdmfAttributeType::GlobalId() );
    nodeIds->setCenter( XdmfAttributeCenter::Node() );
    nodeIds->setName( "NODEID" );
    nodeIds->insert( 0, _nodeIds.data(), _nNodes, 1, 1 );
    shared_ptr< XdmfInformation > nodeIdsInfo = XdmfInformation::New( "ID", "The nodal IDs" );
    nodeIds->insert( nodeIdsInfo );
    grid->insert( nodeIds );

    //Save the nodal coordinates
    shared_ptr< XdmfGeometry > nodeGeometry = XdmfGeometry::New();
    nodeGeometry->setType( XdmfGeometryType::XYZ() );
    nodeGeometry->setName( "Coordinates" );
    nodeGeometry->insert( 0, _coordinates.data(), 3 * _nNodes, 1, 1 );
    shared_ptr< XdmfInformation > nodeGeometryInfo = XdmfInformation::New( "Coordinates", "Coordinates of the nodes in x1, y1, z1, x2, ... format " );
    nodeGeometry->insert( nodeGeometryInfo );
    grid->setGeometry( nodeGeometry );

    //Save the nodesets
    const std::set< BoundaryID > nodesetIds = _mesh_ptr->meshNodesetIds();

    //TODO: Having to loop over the nodes in the nodesets to determine their size seems sub-optimal
    //      but I suspect having to resize the vectors would be worse.
    std::vector< unsigned int > nodesetSizes( nodesetIds.size(), 0 );
    std::vector< std::shared_ptr< std::vector< unsigned int > > > nodesetNodeIds;
    nodesetNodeIds.reserve( nodesetIds.size() );
    std::vector< shared_ptr< XdmfSet > > XDMFnodeSets( nodesetIds.size() );

    unsigned int indx = 0;

    for ( auto it = nodesetIds.begin(); it != nodesetIds.end(); it++ ){

        //Get the size of the nodeset
        for ( auto nit = _mesh_ptr->getMesh().bid_nodes_begin( *it ); nit != _mesh_ptr->getMesh().bid_nodes_end( *it ); nit++ ){
            Node* node = *nit;

            nodesetSizes[ indx ] += 1;
        }

        //Get the ids of the nodes in the nodeset
        nodesetNodeIds.emplace_back( std::shared_ptr< std::vector< unsigned int > >( new std::vector< unsigned int > ) ) ;
        nodesetNodeIds[ indx ]->reserve( nodesetSizes[ indx ] );

        for ( auto nit = _mesh_ptr->getMesh().bid_nodes_begin( *it ); nit != _mesh_ptr->getMesh().bid_nodes_end( *it ); nit++ ){
            Node* node = *nit;

            nodesetNodeIds[ indx ]->push_back( node->id() );
        }

        //Emplace the nodeset IDs into the XDMF domain
        XDMFnodeSets[ indx ] = XdmfSet::New();
        XDMFnodeSets[ indx ]->setType( XdmfSetType::Node() );
        XDMFnodeSets[ indx ]->setName( _mesh_ptr->getBoundaryName( *it ) );
        XDMFnodeSets[ indx ]->insert( 0, nodesetNodeIds[ indx ]->data(), nodesetSizes[ indx ], 1, 1 );
        grid->insert( XDMFnodeSets[ indx ] );

        indx++;

    }

    /*!==================
    | Save the elements |
    ===================*/

    std::vector< unsigned int > _topology;

    unsigned int _nFiniteElements = 0; //We will ignore infinite elements
    unsigned int _topology_size = 0;

    //Iterate through the elements to determine the size of the topology vector
    for ( auto & elem : as_range( local ? _mesh_ptr->activeLocalElementsBegin() : _mesh_ptr->getMesh().active_elements_begin(),
                                  local ? _mesh_ptr->activeLocalElementsEnd() : _mesh_ptr->getMesh().active_elements_end() ) ){

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
        if ( elem->infinite() ){
            continue;
        }
#endif
        _nFiniteElements += 1; //Increment the number of elements

        if( elem->type( ) == HEX8 ){
            _topology_size += 9; //Add space for the topology type and the number of sides
        }
        else{ //Polyhedron in other cases
            _topology_size += 2; //Add space for the topology type and number of sides
    
            for ( unsigned int i = 0; i < elem->n_sides(); i++ ){
                _topology_size += elem->side_ptr(i)->n_nodes() + 1; //Save the number of nodes for a given side
            }
        }
    }

    _elementIds.reserve( _nFiniteElements ); //Reserve the element id vector
    _subdomainIds.reserve( _nFiniteElements ); //Reserve the subdomain id vector
    _topology.reserve( _topology_size ); //Reserve the topology vector

    //Construct the topology vector
    for ( auto & elem : as_range( local ? _mesh_ptr->activeLocalElementsBegin() : _mesh_ptr->getMesh().active_elements_begin(),
                                  local ? _mesh_ptr->activeLocalElementsEnd() : _mesh_ptr->getMesh().active_elements_end() ) ){

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
        if ( elem->infinite() ){
            continue;
        }
#endif
        _elementIds.push_back( elem->id() ); //Save the element ID
        _subdomainIds.push_back( elem->subdomain_id() ); //Save the subdomain ID

        if ( elem->dim() == 1 ){
            _topology.push_back(  2 ); //XDMF Polyline Type
        }
        else if ( elem->dim() == 2 ){
            _topology.push_back(  3 ); //XDMF Polygon Type
        }
        else if ( elem->dim() == 3 ){
            if ( elem->type( ) == HEX8 ){
                _topology.push_back( 9 );
            }
            else{
                _topology.push_back( 16 ); //XDMF Polyhedron Type
            }
        }
        else{
            mooseError( "Element dimension not implemented" );
        }

        if ( elem->type( ) == HEX8 ){

            for ( unsigned int i = 0; i < elem->n_nodes( ); i++ ){

                _topology.push_back( elem->node_id( i ) );

            }

        }
        else{
        _topology.push_back( elem->n_sides() ); //Store the number of sides for the element
    
            for ( unsigned int i = 0; i < elem->n_sides(); i++ ){
    
                if ( elem->dim() > 1 ){
    
                    _topology.push_back( elem->side_ptr( i )->n_nodes() ); //Store the number of nodes on the given side
    
                }
    
                for ( unsigned int n = 0; n < elem->side_ptr( i )->n_nodes(); n++ ){
    
                    _topology.push_back( elem->side_ptr( i )->node_ptr( n )->id() ); //Push back the node id's on this face
    
                }
    
            }
        }
    }

    //Set the topology
    shared_ptr< XdmfTopology > topology = XdmfTopology::New();
    topology->setType( XdmfTopologyType::Mixed() );
    topology->setName( "Topology" );
    topology->insert( 0, _topology.data(), _topology_size, 1, 1 );
    grid->setTopology( topology );

    //Set the element ID numbers
    shared_ptr< XdmfAttribute > elementIds = XdmfAttribute::New( );
    elementIds->setType( XdmfAttributeType::GlobalId( ) );
    elementIds->setCenter( XdmfAttributeCenter::Cell( ) );
    elementIds->setName( "ELEMID" );
    elementIds->insert( 0, _elementIds.data( ), _nFiniteElements, 1, 1 );
    shared_ptr< XdmfInformation > elementIdsInfo = XdmfInformation::New( "ID", "The element IDs" );
    elementIds->insert( elementIdsInfo );
    grid->insert( elementIds );

//    //Set the element subdomain ID numbers
//    shared_ptr< XdmfAttribute > subdomainIds = XdmfAttribute::New( );
//    subdomainIds->setType( XdmfAttributeType::GlobalId( ) );
//    subdomainIds->setCenter( XdmfAttributeCenter::Cell( ) );
//    subdomainIds->setName( "SUBDOMAINID" );
//    subdomainIds->insert( 0, _subdomainIds.data( ), _nFiniteElements, 1, 1 );
//    shared_ptr< XdmfInformation > subdomainIdsInfo = XdmfInformation::New( "BLOCK ID", "The subdomain ( block ) IDs" );
//    subdomainIds->insert( subdomainIdsInfo );
//    grid->insert( subdomainIds );


    //Save the element subdomains

    const std::set< SubdomainID > subdomainIds = _mesh_ptr->meshSubdomains();
    //TODO: Having to loop over the elements to determine the size of the subdomains seems inefficient
    //      If this could be determined without a loop that would be ideal.

    std::vector< unsigned int > subdomainSizes( subdomainIds.size(), 0 );
    std::vector< std::shared_ptr< std::vector< unsigned int > > > subdomainElementIds;
    subdomainElementIds.reserve( subdomainElementIds.size() );
    std::vector< shared_ptr< XdmfSet > > XDMFelementSets( subdomainIds.size() );

    indx = 0;

    for ( auto it = subdomainIds.begin(); it != subdomainIds.end(); it++ ){

//        std::cout << _mesh_ptr->getSubdomainName( *it ) << ": " << *it  << "\n";

        //Get the size of the nodeset
        for ( auto eit =  _mesh_ptr->getMesh().active_subdomain_elements_begin( *it );
                   eit != _mesh_ptr->getMesh().active_subdomain_elements_end( *it ); eit++ ){
            Elem* elem = *eit;

            subdomainSizes[ indx ] += 1;
        }

        //Get the ids of the nodes in the nodeset
        subdomainElementIds.emplace_back( std::shared_ptr< std::vector< unsigned int > >( new std::vector< unsigned int > ) ) ;
        subdomainElementIds[ indx ]->reserve( subdomainSizes[ indx ] );

        for ( auto eit =  _mesh_ptr->getMesh().active_subdomain_elements_begin( *it );
                   eit != _mesh_ptr->getMesh().active_subdomain_elements_end( *it ); eit++ ){
            Elem* elem = *eit;

            subdomainElementIds[ indx ]->push_back( elem->id() );
        }

        //Emplace the nodeset IDs into the XDMF domain
        XDMFelementSets[ indx ] = XdmfSet::New();
        XDMFelementSets[ indx ]->setType( XdmfSetType::Cell() );
        XDMFelementSets[ indx ]->setName( _mesh_ptr->getSubdomainName( *it ) );
        XDMFelementSets[ indx ]->insert( 0, subdomainElementIds[ indx ]->data(), subdomainSizes[ indx ], 1, 1 );
        grid->insert( XDMFelementSets[ indx ] );

        indx++;

    }



    /*==========================
    | Write to the output file |
    ==========================*/

    shared_ptr< XdmfGridCollection > collection = domain->getGridCollection( _num_temporal_collections - 1 );

    collection->insert( grid );

    //Construct the file writers
    shared_ptr< XdmfHDF5Writer > heavyWriter = XdmfHDF5Writer::New( _filename + ".h5" );
    heavyWriter->setReleaseData( true );
    if ( ( _nNodes < 1000 ) || ( _nFiniteElements < 1000 ) ){
        heavyWriter->setChunkSize( std::max( _nNodes, _nFiniteElements ) );
    }
    shared_ptr< XdmfWriter > writer = XdmfWriter::New( _filename + ".xdmf", heavyWriter );
    writer->setLightDataLimit( 1 );

    //Write out the data
    domain->accept( writer );

    //Update the reference grid
    _current_reference_grid = domain->getGridCollection( 0 )->getNumberUnstructuredGrids( ) - 1;

    return;
}

void Xdmf::outputNodalVariables( const std::set< std::string > * system_names ){
    /*!
     * Write the nodal variables out to the XDMF file
     */

    //Collect all of the variable names
    std::vector< std::string > names;
    _es_ptr->build_variable_names( names, nullptr, system_names );

    //Determine the nodal output names
    const std::vector< std::string > outputNames( getNodalVariableOutput().begin(),
                                                  getNodalVariableOutput().end() );

    //Determine the total number of variables and nodes
    unsigned int num_vars = names.size();
    unsigned int num_nodes = _mesh_ptr->getMesh().n_nodes();

    if ( _mesh_ptr->processor_id() ){
        return;
    }

    //Construct the solution vector
    std::vector< Number > soln;
    _es_ptr->build_solution_vector( soln, system_names );

    //Begin looping through the variables
    for ( unsigned int c = 0; c < outputNames.size(); c++ ){

        //Find the current output variable
        std::vector< std::string >::const_iterator pos = std::find( names.begin(), names.end(), outputNames[ c ] );

        if ( pos == names.end() ){
            continue;
        }

        //Determine the current variable index
        unsigned int variable_name_position = cast_int< unsigned int >( pos - names.begin() );

        //Construct the output vectors
#ifdef LIBMESH_USE_REAL_NUMBERS
        std::vector< Number > cur_soln;

        cur_soln.reserve( num_nodes );

#else
        std::vector< Real > real_parts;
        std::vector< Real > imag_parts;
        std::vector< Real > magnitudes;

        real_parts.reserve( num_nodes );
        imag_parts.reserve( num_nodes );
        magnitudes.reserve( num_nodes );
#endif

        //Extract the variable values
        for ( const auto & node : _mesh_ptr->getMesh().node_ptr_range() ){
            
            unsigned int idx = node->id() * num_vars + variable_name_position;
#ifdef LIBMESH_USE_REAL_NUMBERS
            cur_soln.push_back( soln[ idx ] );
#else
            real_parts.push_back( soln[ idx ].real() );
            imag_parts.push_back( soln[ idx ].imag() );
            magnitudes.push_back( std::abs( soln[ idx ] ) );
#endif

        }

        //Write nodal variables to the file
        //Open the output file
        shared_ptr< XdmfReader > reader = XdmfReader::New();
    
        //Read in the XDMF domain    
        shared_ptr< XdmfDomain > domain = shared_dynamic_cast< XdmfDomain >( reader->read( _filename + ".xdmf" ) );
    
        //Get the most recent temporal grid collection
        shared_ptr< XdmfGridCollection > collection = domain->getGridCollection( _num_temporal_collections - 1 );

        //Get the number of grids
        unsigned int nGrids = collection->getNumberUnstructuredGrids( );
    
        if ( nGrids == 0 ){
            mooseError( "No meshes found in XDMF output" );
        }

        //Get the most recent grid //TODO: Should determine how to handle multiple grids
        shared_ptr< XdmfUnstructuredGrid > grid = collection->getUnstructuredGrid( nGrids - 1 );
    
        //Make sure the time associated with the grid is correct
        if ( abs( grid->getTime()->getValue() - _time ) > 1e-9 ){ //TODO: This may fail for very small timesteps
            mooseError( "The time on the most recent grid doesn't coincide with the current time in MOOSE" );
        }

#ifdef LIBMESH_USE_REAL_NUMBERS
        shared_ptr< XdmfAttribute > cur_soln_ptr = XdmfAttribute::New( );
        cur_soln_ptr->setType( XdmfAttributeType::Scalar( ) );
        cur_soln_ptr->setCenter( XdmfAttributeCenter::Node( ) );
        cur_soln_ptr->setName( outputNames[ c ] );

        cur_soln_ptr->insert( 0, cur_soln.data(), num_nodes, 1, 1 );

        shared_ptr< XdmfInformation > cur_soln_info = XdmfInformation::New( outputNames[ c ], "Nodal quantity " + outputNames[ c ] );

        cur_soln_ptr->insert( cur_soln_info );

        grid->insert( cur_soln_ptr );
#else
        shared_ptr< XdmfAttribute > real_ptr = XdmfAttribute::New( );
        shared_ptr< XdmfAttribute > imag_ptr = XdmfAttribute::New( );
        shared_ptr< XdmfAttribute > mag_ptr  = XdmfAttribute::New( );

        real_ptr->setType( XdmfAttributeType::Scalar( ) );
        imag_ptr->setType( XdmfAttributeType::Scalar( ) );
        mag_ptr->setType( XdmfAttributeType::Scalar( ) );

        real_ptr->setCenter( XdmfAttributeCenter::Node( ) );
        imag_ptr->setCenter( XdmfAttributeCenter::Node( ) );
        mag_ptr->setCenter( XdmfAttributeCenter::Node( ) );

        real_ptr->setName( outputNames[ c ] );
        imag_ptr->setName( outputNames[ c ] );
        mag_ptr->setName( outputNames[ c ] );

        real_ptr->insert( 0, real_parts.data(), num_nodes, 1, 1 );
        imag_ptr->insert( 0, imag_parts.data(), num_nodes, 1, 1 );
        mag_ptr->insert( 0, magnitudes.data(), num_nodes, 1, 1 );

        shared_ptr< XdmfInformation > real_info = XdmfInformation::New( "Real " + outputNames[ c ], "Real part of the nodal quantity " + outputNames[ c ] );
        shared_ptr< XdmfInformation > imag_info = XdmfInformation::New( "Imaginary " + outputNames[ c ], "Imaginary part of the nodal quantity " + outputNames[ c ] );
        shared_ptr< XdmfInformation > mag_info = XdmfInformation::New( "Magnitude " + outputNames[ c ], "Magnitude of the nodal quantity " + outputNames[ c ] );

        real_ptr->insert( real_info );
        imag_ptr->insert( imag_info );
        mag_ptr->insert( mag_info );

        grid->insert( real_ptr );
        grid->insert( imag_ptr );
        grid->insert( mag_ptr );
#endif

        //Construct the file writers
        shared_ptr< XdmfHDF5Writer > heavyWriter = XdmfHDF5Writer::New( _filename + ".h5" );
        heavyWriter->setReleaseData( true );
        if ( num_nodes < 1000 ){
            heavyWriter->setChunkSize( num_nodes );
        }
        shared_ptr< XdmfWriter > writer = XdmfWriter::New( _filename + ".xdmf", heavyWriter );
        writer->setLightDataLimit( 1 );
    
        //Write out the data
        domain->accept( writer );
    }

    return;
}

