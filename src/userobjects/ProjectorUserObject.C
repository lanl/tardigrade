//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ProjectorUserObject.h"

registerMooseObject("tardigradeApp", ProjectorUserObject);

template <>
InputParameters
validParams<ProjectorUserObject>()
{
    InputParameters params = validParams<ElementUserObject>();
    params.addRequiredParam<UserObjectName>("nodal_overlap_userobject", "The NodalOverlapUserObject that detects which micro-scale nodes overlap with the macro-scale");
    params.addRequiredParam<UserObjectName>("element_integrate_userobject", "The ElementIntegrateUserObject that computes shape-function weighted integrals at the nodes");

    return params;
}

ProjectorUserObject::ProjectorUserObject(const InputParameters & parameters)
    : ElementUserObject(parameters),
    _nodal_overlap(getUserObject<NodalOverlapUserObject>("nodal_overlap_userobject")),
    _element_integrate(getUserObject<ElementIntegrateUserObject>("element_integrate_userobject"))
{
}

void
ProjectorUserObject::initialize()
{
    _console << "Initializing Projector UserObject: " << name() << std::endl;

    //Resize the shape-function matrix terms to zero
    tripletList.resize(0);

    //Get the macro to row and micro to col maps
    macro_node_to_col = _nodal_overlap.get_macro_node_to_col();
    micro_node_to_row = _nodal_overlap.get_micro_node_to_row();
}

void
ProjectorUserObject::execute()
{

    const std::vector< dof_id_type> *micro_nodes = _nodal_overlap.get_relevant_micro_nodes(_current_elem->id());
    overlap::vecOfvec local_micronode_positions;

    overlap::vecOfvec local_nodes;
    overlap::vecOfvec local_gpts;

    std::vector< dof_id_type > macro_node_ids; //!The id numbers of the nodes associated with the element

    //!Check if any micro-nodes are overlapped by the current element
    if (micro_nodes){
        std::cout << "Found nodes in element\n";

        //Collect the local coordinates of the element's nodes
        collect_local_nodes(local_nodes, macro_node_ids);

        //Collect the locations of the quadrature points
        collect_quadrature_points(local_gpts);

        //Collect the coordinates of the micro-nodes in the element's local coordinates
        collect_local_micronode_positions(local_micronode_positions);

        //Compute the weights for the overlap coupling
        dns_weights.insert( std::pair< dof_id_type, std::vector< overlap::integrateMap > >(_current_elem->id(), std::vector< overlap::integrateMap >(local_gpts.size())));
        compute_dns_weights(local_nodes, local_gpts, micro_nodes, local_micronode_positions);

        //Add contributions to the shapefunction matrix
        compute_shapefunction_matrix_contributions(macro_node_ids);

    }
}

void
ProjectorUserObject::threadJoin(const UserObject & y)
{

}

void
ProjectorUserObject::finalize()
{
    std::cout << "Finalizing Projector UserObject\n";
//    for (unsigned int i=0; i<integrated_weights.size(); i++){
//        std::cout << integrated_weights[i] << ", " << integrated_weighted_densities[i] << ", " << integrated_weighted_densities[i]/integrated_weights[i]  << "\n";
//    }
//    mooseError("derp");
    std::cout << "\tConstructing the shape-function matrix\n";
    shapefunction = overlap::SpMat(n_micro_dof*micro_node_to_row->size(), n_macro_dof*macro_node_to_col->size()); //Size the sparse matrix
    shapefunction.setFromTriplets(tripletList.begin(), tripletList.end()); //Fill the sparse matrix
    shapefunction.makeCompressed(); //Compress the shapefunction matrix

    //Get the micro-information
    unsigned int num_macro_free, num_macro_ghost, num_micro_free, num_micro_ghost;
    _nodal_overlap.get_node_info(num_macro_ghost, num_macro_free, num_micro_ghost, num_micro_free);

    std::cout << "num_macro_free:  " << num_macro_free << "\n";
    std::cout << "num_macro_ghost: " << num_macro_ghost << "\n";
    std::cout << "num_micro_free:  " << num_micro_free << "\n";
    std::cout << "num_micro_ghost: " << num_micro_ghost << "\n";

    //Extract the shape-function matrix sub-blocks
    overlap::SpMat NQD   = shapefunction.block( 0, 0, n_micro_dof*num_micro_free, n_macro_dof*num_macro_free);
    overlap::SpMat NQDh  = shapefunction.block( 0, n_macro_dof*num_macro_free, n_micro_dof*num_micro_free, n_macro_dof*num_macro_ghost);
    overlap::SpMat NQhD  = shapefunction.block( n_micro_dof*num_micro_free, 0, n_micro_dof*num_micro_ghost,  n_macro_dof*num_macro_free);
    overlap::SpMat NQhDh = shapefunction.block( n_micro_dof*num_micro_free, n_macro_dof*num_macro_free, n_micro_dof*num_micro_ghost, n_macro_dof*num_macro_ghost);

//    NQDh.makeCompressed();

    //Compute the Projection Matrices

    //Solve for BDhQ
    std::cout << "Solving for BDhQ\n";
    overlap::SpMat Imicro(n_micro_dof*num_micro_free, n_micro_dof*num_micro_free); //Set the right-hand side
    Imicro.setIdentity();
    solve_for_projector(NQDh, Imicro, BDhQ); //Run the QR solver

    //Solve for BDhD (don't need to do this. Just for generality)
    std::cout << "Solving for BDhD\n";
//    overlap::solve_for_projector(NQDh, NQD, BDhD);
    BDhD = overlap::SpMat(n_macro_dof*num_macro_ghost, n_macro_dof*num_macro_free);

    //Solve for BQhD
    std::cout << "Solving for BQhD\n";
    BQhD = NQhD;// + NQhDh*BDhD;

    //Solve for BQhQ
    std::cout << "Solving for BQhQ\n";
    BQhQ = NQhDh*BDhQ;

//    std::cout << "NQD:\n";
//    std::cout << NQD << "\n\n";
//
//    std::cout << "NQDh:\n";
//    std::cout << NQDh << "\n\n";
//
//    std::cout << "NQhD:\n";
//    std::cout << NQhD << "\n\n";
//
//    std::cout << "NQhDh:\n";
//    std::cout << NQhDh << "\n\n";
//
//    std::cout << "tripletList.size(): " << tripletList.size() << "\n";
//    std::cout << "Exiting Projector UserObject\n\n";
//    assert(1==0);
}

void
ProjectorUserObject::collect_local_nodes(overlap::vecOfvec &local_nodes, std::vector< dof_id_type > &macro_node_ids){
    /*!
    Collect all of the local node coordinates and their id numbers

    :param overlap::vecOfvec local_nodes: The local coordinates of the nodes
    :param std::vector< dof_id_type > macro_node_ids: The id numbers of the nodes
    */

    Point local_node; //Temporary point

    //Reserve the memory required for the local nodes and their ids
    local_nodes.reserve(_current_elem->n_nodes());
    macro_node_ids.reserve(local_nodes.size());

    //Loop through the element's nodes and assign the required values
    for (unsigned int i=0; i<_current_elem->n_nodes(); i++){
        local_node = _current_elem->master_point(i);
        local_nodes.push_back({local_node(0), local_node(1), local_node(2)});
        macro_node_ids.push_back(_current_elem->node_id(i));
    }
}

void
ProjectorUserObject::collect_quadrature_points(overlap::vecOfvec &local_gpts){
    /*!
    Collect the local coordinates of the quadrature points of the element

    :param overlap::vecOfvec local_gpts: The local coordinates of the quadrature points
    */

    //Get the local coordinates of the quadrature points in MOOSE Point format
    const std::vector< Point > local_gpt_positions = _assembly.qRule()->get_points();

    //Reserve the memory for the quadrature points in vector format
    local_gpts.reserve(local_gpt_positions.size());

    //Loop through the quadrature points
    for (unsigned int i=0; i<local_gpt_positions.size(); i++){
        local_gpts.push_back({local_gpt_positions[i](0), local_gpt_positions[i](1), local_gpt_positions[i](2)});
    }

}

void
ProjectorUserObject::collect_local_micronode_positions(overlap::vecOfvec &local_micronode_positions_vec){
    /*!
    Collect the local positions of the DNS in the macro element

    :param overlap::vecOfvec local_micronode_positions_vec: The coordinates of the DNS nodes in the macro-scale element's local coordinates
    */

    //Get the local coordinates of the quadrature points in MOOSE Point format
    const std::vector< Point >* local_micronode_positions = _nodal_overlap.get_local_node_positions( _current_elem->id());

    //Reserve the memory for the local coordinates of the quadrature points in vector format
    local_micronode_positions_vec.reserve(local_micronode_positions->size());

    //Loop through the micro-nodes
    for (unsigned int i=0; i<local_micronode_positions->size(); i++){
        local_micronode_positions_vec.push_back({(*local_micronode_positions)[i](0),
                                                 (*local_micronode_positions)[i](1),
                                                 (*local_micronode_positions)[i](2)});
    }
}

void
ProjectorUserObject::map_integrateMap_values( unsigned int gpt){
    /*!
    Update the dns weights by mapping them from the local to the global coordinate systems

    :param unsigned int gpt: The number of the quadrature point
    */

    overlap::integrateMap::iterator itiM; //!The iterator through the integration map

    dof_id_type elnum = _current_elem->id();

    std::vector< Point > cell_points; //!A vector of points to compute the shape functions and mappings at

    std::vector< double > ones; //!A vector of ones

    unsigned int index; //!An index variable for iterating through the nodes
    std::vector< double > N; //!A normal vector in master coordinates
    double A; //!An area in master coordinates
    overlap::vecOfvec Finv(3); //! The map from global coordinates to local coordinates
    for (unsigned int i=0; i<3; i++){Finv[i].resize(3);}

    cell_points.reserve(3*dns_weights[elnum][gpt].size()); //TODO: Add the number of points to dns_weights rather than this inefficient estimate

    //Extract all of the points to compute the shape functions and mappings at
    for (itiM=dns_weights[elnum][gpt].begin(); itiM != dns_weights[elnum][gpt].end(); itiM++){
        //Add the cell center to the cell_points vector
        cell_points.push_back(Point(itiM->second.coordinates[0],
                                    itiM->second.coordinates[1],
                                    itiM->second.coordinates[2]));

        //Add the centroids of the surface areas to the cell_points vector
        for (unsigned int j=0; j<itiM->second.face_centroids.size(); j++){
            cell_points.push_back(Point(itiM->second.face_centroids[j][0],
                                        itiM->second.face_centroids[j][1],
                                        itiM->second.face_centroids[j][2]));
        }
    }

    //Set the ones vector
    ones = std::vector< double > (cell_points.size(), 1.);

    //Get the finite element assembly
    auto fe = _assembly.getFE(libMesh::FEType(_current_elem->default_order()), _mesh.dimension());

    //Reinit the fe assembly to compute the jacobians and the transformation maps
    fe->reinit(_current_elem, &cell_points, &ones);

    //Get the required transformation values
    std::vector< Point > xyz = fe->get_xyz();
    std::vector< Real > JxW = fe->get_JxW();
    std::vector< Real > dxidx = fe->get_dxidx();
    std::vector< Real > dxidy = fe->get_dxidy();
    std::vector< Real > dxidz = fe->get_dxidz();
    std::vector< Real > detadx = fe->get_detadx();
    std::vector< Real > detady = fe->get_detady();
    std::vector< Real > detadz = fe->get_detadz();
    std::vector< Real > dzetadx = fe->get_dzetadx();
    std::vector< Real > dzetady = fe->get_dzetady();
    std::vector< Real > dzetadz = fe->get_dzetadz();

    //Apply the transformations to construct the integration weights
    index = 0;
    for (itiM=dns_weights[elnum][gpt].begin(); itiM!=dns_weights[elnum][gpt].end(); itiM++){

        //Transform the volume
        itiM->second.volume *= JxW[index];

        //Set the true coordinates
        for (unsigned int j=0; j<3; j++){
            itiM->second.coordinates[j] = xyz[index](j);
        }

        //Compute the volume, density, and center of gravity in the reference configuration
        volumes[_current_elem->id()][gpt] += itiM->second.volume;
        densities[_current_elem->id()][gpt] += _element_integrate.get_nodal_density( itiM->first )*itiM->second.volume;
        for (unsigned int j=0; j<3; j++){
            cgs[_current_elem->id()][gpt][j] += itiM->second.coordinates[j]*_element_integrate.get_nodal_density( itiM->first )*itiM->second.volume;
        }

        //Increment the index
        index++;

        //Transform the normals
        for (unsigned int j=0; j<itiM->second.das.size(); j++){
            //Set the normal
            N = itiM->second.normal(j);
            A = itiM->second.area(j);

            //Set the inverse of the transformation
            Finv[0][0] = dxidx[index];
            Finv[0][1] = dxidy[index];
            Finv[0][2] = dxidz[index];
            Finv[1][0] = detadx[index];
            Finv[1][1] = detady[index];
            Finv[1][2] = detadz[index];
            Finv[2][0] = dzetadx[index];
            Finv[2][1] = dzetady[index];
            Finv[2][2] = dzetadz[index];

            //Apply Nanson's relation to the differential area
            overlap::apply_nansons_relation(N, JxW[index]*A, Finv, itiM->second.das[j]);

            //Set the location of the face centroids
            for (unsigned int k=0; k<3; k++){
                itiM->second.face_centroids[j][k] = xyz[index](k);
            }

            index++;
        }
    }

    //Remove the total volume and mass from the density and cg calculations
    densities[_current_elem->id()][gpt] /= volumes[_current_elem->id()][gpt];
    for (unsigned int j=0; j<3; j++){
        cgs[_current_elem->id()][gpt][j] /= (volumes[_current_elem->id()][gpt]*densities[_current_elem->id()][gpt]);
    }

    //Reinitialize the element to the default settings
    fe->reinit(_current_elem);
}

void
ProjectorUserObject::compute_dns_weights(const overlap::vecOfvec &local_nodes, const overlap::vecOfvec &local_gpts,
                                         const std::vector< dof_id_type >* micro_nodes, const overlap::vecOfvec &local_micronode_positions){
    /*!
    Compute the weights for the integration of a DNS over the encompassing macro-element.

    :param overlap::vecOfvec local_nodes: The local coordinates of the nodes of the macro-element
    :param overlap::vecOfvec local_gpts: The local coordinates of the quadrature points  of the macro-element
    :param std::vector< dof_id_type > micro_nodes: The id numbers of the DNS nodes encompassed by the macro-element
    :param overlap::vecOfvec local_micronode_positions: The local coordinates of the DNS nodes
    */

    dof_id_type elnum = _current_elem->id();

    //Compute the weights for the overlap coupling
    overlap::OverlapCoupling oc = overlap::OverlapCoupling(local_nodes, local_gpts);
    oc.compute_weights(*micro_nodes, local_micronode_positions, dns_weights[elnum]);

    //Add the current element to the map for the volumes, densities, and centers of gravity
    volumes.insert( std::pair< dof_id_type, std::vector< double > >(_current_elem->id(), std::vector< double >(dns_weights[elnum].size(), 0.)));
    densities.insert( std::pair< dof_id_type, std::vector< double > >(_current_elem->id(), std::vector< double >(dns_weights[elnum].size(), 0.)));
    cgs.insert( std::pair< dof_id_type, overlap::vecOfvec >(_current_elem->id(), overlap::vecOfvec(dns_weights[elnum].size(), {0., 0., 0.})));

    //Map the dns weights from the master to the global coordinates
    for (unsigned int i=0; i<dns_weights[elnum].size(); i++){
        //Only perform the map if there are nodes in the quadrature point
        if (dns_weights[elnum][i].size()>0){
            map_integrateMap_values(i);
        }
    }
}

void
ProjectorUserObject::get_cg_local_coordinates(unsigned int gpt, Point &local_cg){
    /*!
    Compute the local coordinates of the center of gravity of the DNS of the current quadrature domain

    :param unsigned int gpt: The current quadrature point number
    :param std::vector< Point > &local_cgs: The local coordinates of the quadrature domain's DNS
    */

    //!Populate the global coordinates of the cg
    Point current_cg(cgs[_current_elem->id()][gpt][0],
                     cgs[_current_elem->id()][gpt][1],
                     cgs[_current_elem->id()][gpt][2]);

    //Compute the local coordinates of the center of gravity
    local_cg = FEInterface::inverse_map( _current_elem->dim(), 
                                         libMesh::FEType(_current_elem->default_order()),
                                         _current_elem,
                                         current_cg);

}

void
ProjectorUserObject::get_cg_phi(unsigned int gpt, overlap::vecOfvec &phis){
    /*!
    Compute the shape function values at the center of gravity

    :param unsigned int gpt: The gauss point number
    :param overlap::vecOfvec phi: The shape function values (nodes, qpts)
    */

    //Get the local coordinates of the center of gravity
    Point local_cg;
    get_cg_local_coordinates(gpt, local_cg);

    //Add the point to the points vector
    std::vector< Point > points;
    points.push_back(local_cg);
    
    //Define the weights vector
    std::vector< double > ones(1, 1.);

    //Get the finite element assembly
    auto fe = _assembly.getFE(libMesh::FEType(_current_elem->default_order()), _mesh.dimension());

    //Re-initialize the fe object
    fe->reinit(_current_elem, &points, &ones);

    //Get the phi values
    phis = fe->get_phi();

    //Re-initialize the element
    fe->reinit(_current_elem);
}

void
ProjectorUserObject::compute_shapefunction_matrix_contributions(const std::vector< dof_id_type > &macro_node_ids){
    /*!
    Compute the contribution of the DNS in the quadrature domains to the shape-function matrix

    :param std::vector< dof_id_type > macro_node_ids: The id numbers of the macro-scale nodes
    */


    dof_id_type elnum = _current_elem->id();

    overlap::vecOfvec phis; //!The shape-function values

//    std::cout << "macro_node_to_col->size(): " << macro_node_to_col->size() << "\n";
//    std::cout << "micro_node_to_row->size(): " << micro_node_to_row->size() << "\n";

    for (unsigned int gpt=0; gpt<dns_weights[elnum].size(); gpt++){
        if (dns_weights[elnum][gpt].size()>0){
            get_cg_phi(gpt, phis);
            overlap::construct_triplet_list(macro_node_to_col, micro_node_to_row, macro_node_ids,
                                            cgs[_current_elem->id()][gpt], phis, dns_weights[elnum][gpt],
                                            tripletList);
        }
    }
}

template< class myType >
void
ProjectorUserObject::print_vector(const std::vector< myType > &v){
    /*!
    Print a vector to the terminal

    :param std::vector< double > &v: The vector to be printed
    */

    for (unsigned int i=0; i<v.size(); i++){
        std::cout << v[i] << " ";
    }
    std::cout << "\n";
}

template< class myType >
void
ProjectorUserObject::print_matrix(const std::vector< std::vector< myType > > &m){
    /*!
    Print a matrix to the terminal

    :param std::vector< std::vector< double > > m: The matrix to be printed
    */

    for (unsigned int i=0; i<m.size(); i++){
        print_vector(m[i]);
    }
}

void
ProjectorUserObject::solve_for_projector(const overlap::SpMat &A, const overlap::SpMat &B, overlap::SpMat &X){
    /*!
    Solve for the projector by solving the sparse matrix linear algebra problem using QR decomposition.

    AX = B

    :param overlap::SpMat A: The left-hand-side matrix
    :param overlap::SpMat B: The right-hand-side matrix
    :param overlap::SpMat X: The solution matrix.
    */

    //Define the solver
    overlap::QRsolver solver;

    //Perform the decomposition
    solver.compute(A);
    if( solver.info() != Eigen::Success){
        std::cout << "Error: Least squares solution to solving for the projector failed\n";
        assert(1==0);
    }
    X = solver.solve(B);
}

const overlap::SpMat* ProjectorUserObject::get_BDhD() const{
    /*!
    Return a pointer to the BDhD matrix
    */

    return &BDhD;
}

const overlap::SpMat* ProjectorUserObject::get_BDhQ() const{
    /*!
    Return a pointer to the BDhQ matrix
    */

    return &BDhQ;
}

const overlap::SpMat* ProjectorUserObject::get_BQhD() const{
    /*!
    Return a pointer to the BQhD matrix
    */

    return &BQhD;
}

const overlap::SpMat* ProjectorUserObject::get_BQhQ() const{
    /*!
    Return a pointer to the BQhQ matrix
    */

    return &BQhQ;
}
