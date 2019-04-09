//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeProjectorAux.h"

registerMooseObject("tardigradeApp", ComputeProjectorAux);

template <>
InputParameters
validParams<ComputeProjectorAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<UserObjectName>("nodal_overlap_userobject", "The NodalOverlapUserObject that detects which micro-scale nodes overlap with the macro-scale");
  params.addRequiredParam<UserObjectName>("element_integrate_userobject", "The ElementIntegrateUserObject that computes shape-function weighted integrals at the nodes");
  return params;
}

ComputeProjectorAux::ComputeProjectorAux(const InputParameters & parameters)
  : AuxKernel(parameters),

  _nodal_overlap(getUserObject<NodalOverlapUserObject>("nodal_overlap_userobject")),
  _element_integrate(getUserObject<ElementIntegrateUserObject>("element_integrate_userobject"))
{
}

/**
 * Auxiliary Kernels override computeValue() instead of computeQpResidual().  Aux Variables
 * are calculated either one per elemenet or one per node depending on whether we declare
 * them as "Elemental (Constant Monomial)" or "Nodal (First Lagrange)".  No changes to the
 * source are necessary to switch from one type or the other.
 */

Real
ComputeProjectorAux::computeValue()
{
    std::ofstream out_file;
    out_file.open("overlap.txt");

    const std::vector< dof_id_type> *micro_nodes = _nodal_overlap.get_relevant_micro_nodes(_current_elem->id());
    const std::vector< Point > * local_node_positions;
    overlap::vecOfvec local_node_positions_vec;
    const std::vector< Point > local_gpt_positions = _assembly.qRule()->get_points();

    overlap::vecOfvec local_nodes;
    overlap::vecOfvec local_gpts;
    std::vector< overlap::integrateMap > dns_weights;
    std::vector< Real > JxW;
    std::vector< Real > dxidx;
    std::vector< Real > dxidy;
    std::vector< Real > dxidz;
    std::vector< Real > detadx;
    std::vector< Real > detady;
    std::vector< Real > detadz;
    std::vector< Real > dzetadx;
    std::vector< Real > dzetady;
    std::vector< Real > dzetadz;
    std::vector< Point > xyz;

    std::vector< double > volumes;
    std::vector< double > densities;
    std::vector< std::vector< double > > cgs;

    unsigned int index;
    std::vector< double > N(3);
    double A;
    overlap::vecOfvec Finv(3);
    for (unsigned int i=0; i<3; i++){Finv[i].resize(3);}
    std::vector< dof_id_type > macro_node_ids;

    //Get the macro to row and micro to col maps
    const std::map< dof_id_type, unsigned int >* macro_node_to_row = _nodal_overlap.get_macro_node_to_row();
    const std::map< dof_id_type, unsigned int >* micro_node_to_col = _nodal_overlap.get_micro_node_to_col();

    //Initialize the sparse matrix
    std::cout << "initializing the sparse matrix...\n";
    overlap::SpMat shapefunction(12*macro_node_to_row->size(), 3*micro_node_to_col->size());
    overlap::SpMat sub_shapefunction(12*macro_node_to_row->size(), 3*micro_node_to_col->size());
    std::cout << "we did it!\n";

    //!Get the micro-nodes from the nodal overlap object
    if (micro_nodes){
        std::cout << "Found nodes in element\n";

        Point local_node;
        Point global_node;

        local_nodes.reserve(_current_elem->n_nodes());
        macro_node_ids.reserve(local_nodes.size());
        for (unsigned int i=0; i<_current_elem->n_nodes(); i++){
            local_node = _current_elem->master_point(i);

            local_nodes.push_back({local_node(0), local_node(1), local_node(2)});
            
            macro_node_ids.push_back(_current_elem->node_id(i));
        }


        //!Output the quadrature point positions
        local_gpts.reserve(local_gpt_positions.size());
        for (unsigned int i=0; i<local_gpt_positions.size(); i++){
            local_gpts.push_back({local_gpt_positions[i](0), local_gpt_positions[i](1), local_gpt_positions[i](2)});
        }

        std::cout << "Assemble micro node local position vector\n";
        local_node_positions = _nodal_overlap.get_local_node_positions( _current_elem->id());

        local_node_positions_vec.reserve(local_node_positions->size());
        for (unsigned int i=0; i<local_node_positions->size(); i++){
            local_node_positions_vec.push_back({(*local_node_positions)[i](0),
                                                (*local_node_positions)[i](1),
                                                (*local_node_positions)[i](2)});
        }

        //Get the finite element assembly
        auto fe = _assembly.getFE(libMesh::FEType(_current_elem->default_order()), _mesh.dimension());

        //Compute the weights for the overlap coupling
        overlap::OverlapCoupling oc = overlap::OverlapCoupling(local_nodes, local_gpts);
        oc.compute_weights(*micro_nodes, local_node_positions_vec, dns_weights);
        overlap::integrateMap::iterator itiM;
        std::vector< Point > cell_points;
        std::vector< double > ones;

        volumes = std::vector< double >(dns_weights.size(), 0.);
        densities = std::vector< double >(dns_weights.size(), 0.);
        cgs = std::vector< std::vector< double > >(dns_weights.size(), {0., 0., 0.});

        for (unsigned int i=0; i<dns_weights.size(); i++){
            std::cout << "dns_weights[" << i << "].size(): " << dns_weights[i].size() << "\n";
            //Extract all of the locations of the points required
            if (dns_weights[i].size()>0){
                cell_points.resize(0); //TODO: This is probably not as efficient as could be desired
                cell_points.reserve(3*dns_weights[i].size());

                //Get all of the points we require the shape functions and mappings at
                for (itiM=dns_weights[i].begin(); itiM!=dns_weights[i].end(); itiM++){
                    //Add the cell center to the cell points object.
                    cell_points.push_back(Point(itiM->second.coordinates[0],
                                                itiM->second.coordinates[1],
                                                itiM->second.coordinates[2]));
                    //Add the centroids of the surface areas to the cell points object
                    for (unsigned int j=0; j<itiM->second.face_centroids.size(); j++){
                        cell_points.push_back(Point(itiM->second.face_centroids[j][0],
                                                    itiM->second.face_centroids[j][1],
                                                    itiM->second.face_centroids[j][2]));
                    }
                }

                //Set the ones vector
                ones = std::vector< double >(cell_points.size(), 1.);

                //Reinit the fe assembly to compute the jacobians and transformation maps
                fe->reinit(_current_elem, &cell_points, &ones);

                //Get the required values
                xyz = fe->get_xyz();
                JxW = fe->get_JxW();
                dxidx = fe->get_dxidx();
                dxidy = fe->get_dxidy();
                dxidz = fe->get_dxidz();
                detadx = fe->get_detadx();
                detady = fe->get_detady();
                detadz = fe->get_detadz();
                dzetadx = fe->get_dzetadx();
                dzetady = fe->get_dzetady();
                dzetadz = fe->get_dzetadz();

                //Apply the transformations to construct the integration weights
                index = 0;
                for (itiM=dns_weights[i].begin(); itiM!=dns_weights[i].end(); itiM++){

                    //Transform the volume
                    itiM->second.volume *= JxW[index];

                    //Set the true coordinates
                    for (unsigned int j=0; j<3; j++){
                        itiM->second.coordinates[j] = xyz[index](j);
                    }

                    //Compute the volume, density, and center of gravity in the reference configuration
                    volumes[i] += itiM->second.volume;
                    densities[i] += _element_integrate.get_nodal_density( itiM->first )*itiM->second.volume;
                    for (unsigned int j=0; j<3; j++){
                        cgs[i][j] += itiM->second.coordinates[j]*_element_integrate.get_nodal_density( itiM->first )*itiM->second.volume;
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

                //Remove the total volume and mass from the density and cg calcuations
                densities[i] /= volumes[i];
                for (unsigned int j=0; j<3; j++){
                    cgs[i][j] /= (volumes[i]*densities[i]);
                }

                //Get the local coordinates of the cg
                std::vector< Point > local_cg;
                Point current_cg(cgs[i][0], cgs[i][1], cgs[i][2]);
                local_cg.push_back(FEInterface::inverse_map( _current_elem->dim(),
                                                             libMesh::FEType(_current_elem->default_order()),
                                                             _current_elem,
                                                             current_cg));
                std::cout << "local_cg: " << local_cg[0](0) << " " << local_cg[0](1) << " " << local_cg[0](2) << "\n";

                //Reset the size of ones
                ones = std::vector< double >(1,1.);

                //Get the value of the shape function at the center of gravity
                fe->reinit(_current_elem, &local_cg, &ones);
                std::vector< std::vector< double > > phis = fe->get_phi();

                std::cout << "shape functions at cg:\n";
                for (unsigned int j=0; j<phis.size(); j++){
                    std::cout << "  node " << j << ": " << phis[j][0] << "\n";
                }

                std::vector< overlap::T > tripletList;

                //construct the triplet list to populate the sub shape-function matrix
                overlap::construct_triplet_list(macro_node_to_row, micro_node_to_col, macro_node_ids,
                                                cgs[i], phis, dns_weights[i], tripletList);
                std::cout << "tripletList.size(): " << tripletList.size() << "\n";
                sub_shapefunction.setFromTriplets(tripletList.begin(), tripletList.end());

                //Add the values to the full shape-function matrix
                shapefunction += sub_shapefunction;
            }
        }

        //Reset the element
        fe->reinit(_current_elem);

        out_file.close();

        mooseError("derp");

    }


    
    
    return 0;
}
