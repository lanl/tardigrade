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
    std::vector< double > nodal_density(0);
    std::vector< double > nodal_volume(0);

    overlap::vecOfvec local_nodes;
    overlap::vecOfvec local_gpts;
    std::vector< overlap::integrateMap > dns_local_weights;
    std::vector< Real > JxW;
    std::vector< Real > dxidx;
    std::vector< Real > dxidy;
    std::vector< Real > dxixz;
    std::vector< Real > detadx;
    std::vector< Real > detady;
    std::vector< Real > detaxz;
    std::vector< Real > dzetadx;
    std::vector< Real > dzetady;
    std::vector< Real > dzetaxz;

    std::vector< double > N(3);
    overlap::vecOfvec Finv(3);
    for (unsigned int i=0; i<3; i++){Finv[i].resize(3);}

    //!Get the micro-nodes from the nodal overlap object
    if (micro_nodes){
        std::cout << "Found nodes in element\n";

        Point local_node;
        Point global_node;

        local_nodes.reserve(_current_elem->n_nodes());
        for (unsigned int i=0; i<_current_elem->n_nodes(); i++){
            local_node = _current_elem->master_point(i);

            local_nodes.push_back({local_node(0), local_node(1), local_node(2)});
            
        }


        //!Output the quadrature point positions
        local_gpts.reserve(local_gpt_positions.size());
        for (unsigned int i=0; i<local_gpt_positions.size(); i++){
//            out_file << local_gpt_positions[i](0) << " " << local_gpt_positions[i](1) << " " << local_gpt_positions[i](2) << "\n";
            local_gpts.push_back({local_gpt_positions[i](0), local_gpt_positions[i](1), local_gpt_positions[i](2)});
        }

        std::cout << "Assemble node local position vector\n";
        nodal_density.reserve(micro_nodes->size());
        nodal_volume.reserve(micro_nodes->size());
        local_node_positions = _nodal_overlap.get_local_node_positions( _current_elem->id());

        local_node_positions_vec.reserve(local_node_positions->size());
        for (unsigned int i=0; i<local_node_positions->size(); i++){
            local_node_positions_vec.push_back({(*local_node_positions)[i](0),
                                                (*local_node_positions)[i](1),
                                                (*local_node_positions)[i](2)});
        }

        const std::vector< Real > local_weights(local_node_positions->size(), 1);
        const std::vector< Real > *local_weights_ptr;
        local_weights_ptr = &local_weights;

        for (unsigned int i=0; i<micro_nodes->size(); i++){
            nodal_density.push_back(_element_integrate.get_nodal_density( (*micro_nodes)[i]));
            nodal_volume.push_back(_element_integrate.get_nodal_volume( (*micro_nodes)[i]));
        }

        //Get the finite element assembly
        auto fe = _assembly.getFE(libMesh::FEType(_current_elem->default_order()), _mesh.dimension());

        //Compute the weights for the overlap coupling
        overlap::OverlapCoupling oc = overlap::OverlapCoupling(local_nodes, local_gpts);
        oc.compute_weights(*micro_nodes, local_node_positions_vec, dns_local_weights);
        overlap::integrateMap::iterator itiM;
        std::vector< Point > cell_points;
        std::vector< double > ones;
        unsigned int indx=0;
        for (unsigned int i=0; i<dns_local_weights.size(); i++){

            //Extract all of the locations of the points required
            cell_points.resize(dns_local_weights[i].size());
            for (itiM=dns_local_weights[i].begin(); itiM!=dns_local_weights[i].end(); itiM++){
                //Add the cell center to the cell points object.
                cell_points.push_back(Point(itiM->second.coordinates[0],
                                            itiM->second.coordinates[1],
                                            itiM->second.coordinates[2]));
                //Add the centroids of the surface areas to the cell points object
                for (unsigned int j=0; j<itiM->second.areas.size(); j++){
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

            //Apply the transformations
            index = 0;
            for (itiM-dns_local_weights[i].begin(); itiM!=dns_local_weights[i].end(); itiM++){
                //Transform the volume
                dns_local_weights[i][itiM->first].volume *= JxW[index]; index++;
                //Transform the normals
                for (unsigned int j=0; j<itiM->second.normals.size(); j++){
                    //Set the normal
                    N = itiM->second.normals[j];

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

                    //Set the normals to the new values and scale them by the area and jacobian
                    itiM->second.normals[j][0] = JxW[index]*itiM->second.areas[j]*(Finv[0][0]*N[0] + Finv[1][0]*N[1] + Finv[2][0]*N[2]);
                    itiM->second.normals[j][1] = JxW[index]*itiM->second.areas[j]*(Finv[0][1]*N[0] + Finv[1][1]*N[1] + Finv[2][1]*N[2]);
                    itiM->second.normals[j][2] = JxW[index]*itiM->second.areas[j]*(Finv[0][2]*N[0] + Finv[1][2]*N[1] + Finv[2][2]*N[2]);
                    index++;
                }
            }
        }

        //Reset the element
        fe->reinit(_current_elem);

        std::cout << "local_node_positions->size(): " << local_node_positions->size() << "\n";
        std::cout << "Micro-nodes\n";


        for (unsigned int i=0; i<local_node_positions->size(); i++){
//            _element_integrate.get_nodal_density( micro_nodes[i], density);
//            _element_integrate.get_nodal_volume( micro_nodes[i], volume);

            out_file << (*micro_nodes)[i];

            out_file << " " << nodal_volume[i];

            out_file << " " << nodal_density[i];

            out_file << " " << (*local_node_positions)[i](0) << " " << (*local_node_positions)[i](1) << " " << (*local_node_positions)[i](2);

            out_file << "\n";
        }

        out_file.close();
//        fe->reinit(_current_elem);
//        JxW = fe->get_JxW();
//        for (unsigned int i=0; i<JxW.size(); i++){
//            std::cout << "JxW[" << i << "]: " << JxW[i] << "\n";
//        }
        mooseError("derp");

    }


    
    
    return 0;
}
