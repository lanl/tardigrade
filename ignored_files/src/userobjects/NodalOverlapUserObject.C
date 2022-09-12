//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NodalOverlapUserObject.h"

//#include "libmesh/bounding_box.h"
//#include "libmesh/mesh_tools.h"

registerMooseObject("tardigradeApp", NodalOverlapUserObject);

template <>
InputParameters
validParams<NodalOverlapUserObject>()
{
    InputParameters params = validParams<NodalUserObject>();
    params.addRequiredParam<SubdomainName>("macroscale_domain", "The micromorphic macroscale domain");
    return params;
}

NodalOverlapUserObject::NodalOverlapUserObject(const InputParameters & parameters)
    : NodalUserObject(parameters),
    _macroscale_domain(getParam<SubdomainName>("macroscale_domain"))
  
{
}

void
NodalOverlapUserObject::initialize()
{
    mooseWarning("Initializing NodalOverlap");
    _console << "Initializing Nodal Overlap UserObject: " << name() << std::endl;

    //!Set the id of the macro-scale domain
    _macro_id = _mesh.getSubdomainID(_macroscale_domain);
    
//    //!Set the iterators for the elements in the macro-scale domain
//    _macro_element_begin = _mesh.getMesh().active_subdomain_elements_begin(_macro_id);
//    _macro_element_end   = _mesh.getMesh().active_subdomain_elements_end(_macro_id);

    //!Create a bounding box that encompasses the macro-scale domain
    _macro_bounding_box = MeshTools::create_subdomain_bounding_box(_mesh, _macro_id);

    //!Output the bounds to the terminal
    std::cout << "\tMacro-scale bounding box:\n";
    std::cout << "\t\tminimum: ";
    _macro_bounding_box.min().write_unformatted(std::cout);
    std::cout << "\t\tmaximum: ";
    _macro_bounding_box.max().write_unformatted(std::cout);
}

void
NodalOverlapUserObject::execute()
{
    //!Check if the current node is located within the macro-scale's bounding box
    const Point *_current_point = static_cast< const Point* >(_current_node);

    //!Form the map from the micro-scale nodes to their elements. Note: Forming it will only be done once.
    const std::map< dof_id_type, std::vector< dof_id_type > > & micro_node_to_micro_elements = _mesh.nodeToElemMap(); //!A map from the micro-scale node id to the micro-scale element ids

    unsigned int num_elem = 0;

    if (_macro_bounding_box.contains_point(*_current_point))
    {
        //!Since the current node is in the bounding box, we can check if the node is in any of the macro-scale elements
        for (const auto &_macro_element : 
             as_range(_mesh.getMesh().active_subdomain_elements_begin(_macro_id),
                     _mesh.getMesh().active_subdomain_elements_end(_macro_id))){

            //!Check if the element contains the node
            if (_macro_element->contains_point(*_current_point)){
               
                //!Compute and store the local coordinates of the current node 
                Point local_point = FEInterface::inverse_map(_macro_element->dim(),
                                                             libMesh::FEType(_macro_element->default_order()),
                                                             _macro_element,
                                                             *_current_point);
//                std::vector< double > local_coords(LIBMESH_DIM, 0);
//                for (unsigned int i=0; i<local_coords.size(); i++){local_coords[i] = local_point(i);}
//                NodeLocalCoordinates.insert( std::pair< dof_id_type, Point >(_current_node->id(), local_point));

                //!Associate the micro-node with a macro-scale element
                updateMacroMicroMap(_macro_element->id(), _current_node->id(), local_point);

                //!Associate the macro-element with its neighbors
                updateMacroNeighborMap(_macro_element);

                //!Store the macro-scale node ids associated with the macro-element
                updateMacroNodeIds(_macro_element);

                //!Store the node in the id to shape-function col index map
                updateMicroNodeRowMap(_current_node->id());                

                //!Store the macro-element's nodes in the id to shape-function row map
                updateMacroNodeColMap(_macro_element);

                //!Find the micro-elements associated with a given micro-node
                auto node_elements = micro_node_to_micro_elements.find(_current_node->id());
                for (auto elem_id = node_elements->second.begin(); 
                     elem_id != node_elements->second.end();
                     elem_id++){

                    //!Update micro_elements if required
                    updateMicroElements(_macro_element->id(), *elem_id);

                }
                num_elem++; //Increment the number of elements the node is shared between
                if (restrict_nodes_to_one_element){
                    break; //A node can only appear in a single element
                }
            }
        }

        //Set the number of elements the micro-node is contained within
        if (num_elem>1){
            micro_node_elcount.insert( std::pair< dof_id_type, unsigned int >(_current_node->id(), num_elem));
        }
    }
}

void
NodalOverlapUserObject::threadJoin(const UserObject & y)
{
    const NodalOverlapUserObject & pps = static_cast<const NodalOverlapUserObject &>(y);

    std::map< dof_id_type, std::vector< dof_id_type > >::iterator it;

    //Combine the macro_element_to_micro_nodes maps
    for (auto const& _macro_element : pps.macro_element_to_micro_nodes){

        it = macro_element_to_micro_nodes.find(_macro_element.first);
        if (it != macro_element_to_micro_nodes.end()){
            it->second.insert(it->second.end(), _macro_element.second.begin(), _macro_element.second.end());
        }
        else{
            macro_element_to_micro_nodes.insert( std::pair< dof_id_type, std::vector< dof_id_type > >(_macro_element.first,
                                                                                                      _macro_element.second));
        }
    }

    //Combine the micro_element vectors
    for (auto const& micro_elem : pps.micro_elements){

        auto elem_location = micro_elements.find(micro_elem.first);
        if (elem_location == micro_elements.end()){
            micro_elements.insert(std::pair< dof_id_type, std::vector< dof_id_type > >(micro_elem.first, micro_elem.second));
        }
        else{
            elem_location->second.insert( elem_location->second.end(), micro_elem.second.begin(), micro_elem.second.end());
        }
    }
}

void
NodalOverlapUserObject::finalize()
{

    std::cout << "Finalizing Nodal Overlap\n";
    std::cout << "\tNumber of macro-scale elements overlapping microscale domain: " << macro_element_to_micro_nodes.size() << "\n";

    std::map< dof_id_type, std::vector< dof_id_type > >::iterator it1;
    std::map< dof_id_type, unsigned int>::iterator it2;

//    unsigned int num_macro_ghost=0;
//    unsigned int num_micro_free=0;
    unsigned int total_macro_nodes = macro_node_to_col.size();
    unsigned int total_micro_nodes = micro_node_to_row.size();
//    unsigned int num_micro_ghost;

    for (it1=macro_element_to_micro_nodes.begin(); it1 != macro_element_to_micro_nodes.end(); it1++){
        std::cout << "\t  Macro element id: " << it1->first << "\n";
        std::cout << "\t    Number of nodes: " << it1->second.size() << "\n";
//        std::cout << "\t    Nodes: "; print_vector(macro_element_nodes[it1->first]);
        std::cout << "\t    Neighbors: "; print_vector(macro_element_neighbors[it1->first]);
        //!Check the macro-element depths i.e. how many recursive sets of their neighbors deep into the DNS are they
        checkElementDepth(it1->first, element_depth);
        std::cout << "\t    Element Depth: " << element_depth[it1->first] << "\n";

        //Add to the total number of ghost macro-nodes and free micro-nodes
        if (element_depth[it1->first] >= ghost_depth){
            num_macro_ghost += macro_element_nodes[it1->first].size();//Define the macro-nodes contained within this element as ghost
            num_micro_free += it1->second.size(); //Define the micro-nodes contained within this element as free
        }
    }
    num_macro_free = total_macro_nodes - num_macro_ghost;
    num_micro_ghost = total_micro_nodes - num_micro_free;

    std::cout << "Nodal status (counts):\n";
    std::cout << "  Macro (micromorphic) free nodes:  " << num_macro_free << "\n";
    std::cout << "  Macro (micromorphic) ghost nodes: " << num_macro_ghost << "\n";
    std::cout << "  Micro (DNS) free nodes:           " << num_micro_free << "\n";
    std::cout << "  Micro (DNS) ghost nodes:          " << num_micro_ghost << "\n";

    std::cout << "\tNumber of micro-elements associated with overlap:             " << micro_elements.size() << "\n";

//    std::cout << "num_macro_ghost: " << num_macro_ghost << "\n";
//    std::cout << "num_micro_free: " << num_micro_free << "\n";

    //!Define the ordering of the macro-node to the row in the shape-function matrix
    //!The matrix is ordered: | (macro free, micro free)  : (macro free, micro ghost)  |
    //!                       ----------------------------------------------------------
    //!                       | (macro ghost, micro free) : (macro ghost, micro ghost) |

    //Assign the ghost macro-nodes and free-micro nodes
    unsigned int ghost_macro_index = num_macro_free;
    unsigned int free_micro_index = 0;

    unsigned int micro_index_true_zero;

    for (it1 = macro_element_to_micro_nodes.begin(); it1 != macro_element_to_micro_nodes.end(); it1++){

        if (element_depth[it1->first] >= ghost_depth){
        
            //!Assign the macro-ghost index
            for (unsigned int n=0; n<macro_element_nodes[it1->first].size(); n++){
                macro_node_to_col[macro_element_nodes[it1->first][n]] = ghost_macro_index;
                ghost_macro_index++;
            }

            //!Assign the micro-free index
            for (unsigned int n=0; n<it1->second.size(); n++){
                micro_node_to_row[it1->second[n]] = free_micro_index;
                if (free_micro_index == 0){micro_index_true_zero = it1->second[n];} //Store the true first index of the matrix
                free_micro_index++;
            }
        }
    }

//    std::cout << "micro_index_true_zero: " << micro_index_true_zero << "\n";

    //Assign the free macro-nodes and ghost-micro nodes
    unsigned int free_macro_index = 0;
    unsigned int ghost_micro_index = num_micro_free;

    for (it1 = macro_element_to_micro_nodes.begin(); it1 != macro_element_to_micro_nodes.end(); it1++){
        
        if (element_depth[it1->first] < ghost_depth){

            //!Assign the macro-free index
            for (unsigned int n=0; n<macro_element_nodes[it1->first].size(); n++){
                if (macro_node_to_col[macro_element_nodes[it1->first][n]] == 0){
                    macro_node_to_col[macro_element_nodes[it1->first][n]] = free_macro_index;
                    free_macro_index++;
                }
            }

            //!Assign the micro ghost index
            for (unsigned int n=0; n<it1->second.size(); n++){
                if ((micro_node_to_row[it1->second[n]] == 0) && (it1->second[n] != micro_index_true_zero)){
                    micro_node_to_row[it1->second[n]] = ghost_micro_index;
                    ghost_micro_index++;
                }

                //Decrement the number of elements a micro-node is contained within if it is shared between a free and ghost macro-element
                if ((!share_ghost_boundary_nodes) && //Check if nodes shared between ghost and free macro-boundaries should be shared
                    (micro_node_elcount.find(it1->second[n]) != micro_node_elcount.end()) && //Check if the micro-node is on a boundary
                    (micro_node_to_row[it1->second[n]] < num_micro_free)) { //Check if the micro-node is free
                    micro_node_elcount[it1->second[n]] -= 1; //Decrement the number of shared elements
                    if (micro_node_elcount[it1->second[n]] <= 1){
                        micro_node_elcount.erase(it1->second[n]); //Delete any nodes from the map which are only on the border of a single macro-element
                    }
                }
            }
        }
    }

    //Data check
    bool run_data_check = false;
    if (run_data_check){
        //Check that all indices are accounted for
        bool found_index = false;
        for (unsigned int i=0; i<micro_node_to_row.size(); i++){
            found_index = false;
            for (it2 = micro_node_to_row.begin(); it2 != micro_node_to_row.end(); it2++){
                if (it2->second == i){
                    found_index = true;
                    break;
                }
            }
            if (!found_index){
                mooseError("Error: index not found in micro_node_to_row map");
            }
        }

        for (unsigned int i=0; i<macro_node_to_col.size(); i++){
            found_index = false;
            for (it2 = macro_node_to_col.begin(); it2 != macro_node_to_col.end(); it2++){
                if (it2->second == i){
                    found_index = true;
                    break;
                }
            }
            if (!found_index){
                mooseError("Error: index not found in macro_node_to_cow map");
            }
        }

        //Check that all nodes in the ghost micromorphic elements are ghost
        for (it1 = macro_element_nodes.begin(); it1 != macro_element_nodes.end(); it1++){
            std::cout << "element depth: " << element_depth[it1->first] << "\n";
            if (element_depth[it1->first] >= ghost_depth){
                std::cout << "Macro-element " << it1->first << " has " << it1->second.size() << " ghost nodes.\n";
                for (unsigned int n=0; n<it1->second.size(); n++){
                    if (macro_node_to_col[it1->second[n]] < num_macro_free){
                        mooseError("Error: node in ghost macro nodes is free");
                    }
                }
            }
        }

        //Check that all of the nodes in the free DNS elements are free
        for (it1= macro_element_to_micro_nodes.begin(); it1 != macro_element_to_micro_nodes.end(); it1++){
            std::cout << "element depth: " << element_depth[it1->first] << "\n";
            if (element_depth[it1->first] >= ghost_depth){
                std::cout << "Macro-element " << it1->first << " contains " << it1->second.size() << " free DNS nodes.\n";
                for (unsigned int n=0; n<it1->second.size(); n++){
                    if (micro_node_to_row[it1->second[n]] >= num_micro_free){
                        mooseError("Error: node in free DNS nodes is ghost");
                    }
                }
            }
        }

        mooseError("Ending datacheck: All tests successful");
    }

//    std::cout << "micro_node_to_row:\n";
//    for (it2 = micro_node_to_row.begin(); it2!=micro_node_to_row.end(); it2++){
//        std::cout << it2->first << ", " << it2->second << "\n";
//    }
//
//    std::cout << "macro_node_to_col:\n";
//    for (it2 = macro_node_to_col.begin(); it2!=macro_node_to_col.end(); it2++){
//        std::cout << it2->first << ", " << it2->second << "\n";
//    }
//    mooseError("Narf!");

    //Check that the indices are sensible to make sure nothing has gone horribly wrong
    if (ghost_macro_index != total_macro_nodes){
        std::cout << "Error: The final macro ghost index is not equal to the number of macro nodes\n";
        assert(1==0);
    }

    if (free_macro_index != num_macro_free){
        std::cout << "Error: The final macro free index is not equal to the number of macro free nodes\n";
        assert(1==0);
    }

    if (ghost_micro_index != total_micro_nodes){
        std::cout << "Error: The final micro ghost index is not equal to the number of micro nodes\n";
        assert(1==0);
    }

    if (free_micro_index != num_micro_free){
        std::cout << "Error: The final micro free index is not equal to the number of micro free nodes\n";
        assert(1==0);
    }

    std::cout << "End of Nodal Overlap\n\n";
}

void
NodalOverlapUserObject::updateMacroMicroMap(dof_id_type macro_elem_id, dof_id_type node_id, Point local_position)
{
    /*!
    Add an element-node pair to the macro_element_to_micro_node map
    */

    std::map< dof_id_type, std::vector< dof_id_type > >::iterator it = macro_element_to_micro_nodes.find(macro_elem_id);
    if (it != macro_element_to_micro_nodes.end()){
        it->second.push_back(_current_node->id());
        macro_element_to_micro_positions[it->first].push_back(local_position);
    }
    else{
        macro_element_to_micro_nodes.insert( std::pair< dof_id_type, std::vector< dof_id_type > >(macro_elem_id, {node_id}));
        macro_element_to_micro_positions.insert( std::pair< dof_id_type, std::vector< Point > >(macro_elem_id, {local_position}));
    }
}

void
NodalOverlapUserObject::updateMacroNeighborMap(const Elem* macro_elem){
    /*!
    Update the macro-element to neighbors map

    :param Elem* macro_elem: The pointer to the macro element in question.
    */

    const Elem* neighbor;
    std::vector< dof_id_type > neighbors;
    std::map< dof_id_type, std::vector< dof_id_type > >::const_iterator it = macro_element_neighbors.find(macro_elem->id());
    if (it == macro_element_neighbors.end()){
//        std::cout << "new element\n";
        neighbors.reserve(macro_elem->n_neighbors());
        for (unsigned int i=0; i<macro_elem->n_neighbors(); i++){
//            std::cout << "i: " << i << "\n";
            neighbor = macro_elem->neighbor_ptr(i);
            if (neighbor != NULL){
                neighbors.push_back(neighbor->id());
            }
        }
//        std::cout << "neighbors: "; print_vector(neighbors);
        macro_element_neighbors.insert( std::pair< dof_id_type, std::vector< dof_id_type > >(macro_elem->id(), neighbors));
    }
}

void
NodalOverlapUserObject::updateMicroElements(dof_id_type macro_elem_id, dof_id_type micro_elem_id)
{
    /*!
    Update micro_elements with the provided elem_id if it is new.
    */

    //!Search for the element id 
    std::map< dof_id_type, std::vector< dof_id_type> >::iterator it = micro_elements.find(micro_elem_id);
    
    //!Store the element id if it hasn't been identified yet
    if (it == micro_elements.end()){
        micro_elements.insert(std::pair< dof_id_type, std::vector< dof_id_type > >(micro_elem_id, {macro_elem_id}));
    }
    else{
        it->second.push_back(macro_elem_id);
    }
}

const std::vector< dof_id_type >*
NodalOverlapUserObject::get_relevant_micro_nodes(dof_id_type macro_elem_id) const
{
    /*!

    Return the macro elements which contain nodes associated with the micro element

    :param dof_id_type micro_elem_id: The id of the micro element
    :param std::vector< dof_id_type >: The returned vector of macro-element associated with the micro-element

    */

    auto it = macro_element_to_micro_nodes.find(macro_elem_id);
    if (it == macro_element_to_micro_nodes.end()){
        return NULL;
    }
    else{
        return &it->second;
    }

    return NULL;
}

const std::vector< Point >*
NodalOverlapUserObject::get_local_node_positions(dof_id_type macro_element_id) const
{
    /*!
    Return the local positions of the nodes contained by the indicated macro element

    :param dof_id_type macro_element_id: The id number of the macro element
    */

    auto it = macro_element_to_micro_positions.find(macro_element_id);
    if (it == macro_element_to_micro_positions.end()){
        return NULL;
    }
    else{
        return &it->second;
    }
    return NULL;
}

void NodalOverlapUserObject::updateMicroNodeRowMap(const dof_id_type micro_node_id){
    /*!
    Add the indicated node id to the id to shape-function matrix row map. The row number will be 
    determined at finalize. Note that the number is the order. There will be a scale factor applied due to 
    the number of degrees of freedom to be mapped from the micro to macro scales.


    :param dof_id_type micro_node_id: The node id to add
    */

//    std::map< dof_id_type, unsigned int >::iterator it;
//    it = micro_node_to_row.find(micro_node_id);
//    if (it == micro_node_to_row.end()){
//        micro_node_to_row.insert( std::pair< dof_id_type, unsigned int>(micro_node_id, 0));
//    }

    micro_node_to_row.insert( std::pair< dof_id_type, unsigned int>(micro_node_id, 0));
    return;
}

void NodalOverlapUserObject::updateMacroNodeIds(const Elem* macro_element){
    /*!
    Add a count of the number of nodes associated with a given element to the node count map

    :param Elem* macro_element: The pointer to the macro-scale element
    */

    std::vector< dof_id_type > nodes;

    //Add the number of nodes for this macro element to the map if required
    std::map< dof_id_type, std::vector< unsigned int> >::iterator it = macro_element_nodes.find(macro_element->id());
    if (it == macro_element_nodes.end()){
        
        nodes.reserve(macro_element->n_nodes());

        for (unsigned int n=0; n<macro_element->n_nodes(); n++){
            nodes.push_back(macro_element->node_id(n));
        }
        macro_element_nodes.insert( std::pair< dof_id_type, std::vector< unsigned int > >(macro_element->id(), nodes));
    }
}

void NodalOverlapUserObject::updateMacroNodeColMap(const Elem* macro_element){
    /*!
    Add the nodes of the macro element to the id to shape-function matrix column map. The column number will 
    be determined at finalize. Note that the number is the order. There will be a scale factor applied due to 
    the number of degrees of freedom to be mapped from the micro to macro scales.

    :param Elem* macro_element: The pointer to the macro-scale element.
    */

    std::map< dof_id_type, unsigned int >::iterator it;


    for (unsigned int n=0; n<macro_element->n_nodes(); n++){
        it = macro_node_to_col.find(macro_element->node_id(n));
        if (it == macro_node_to_col.end()){
            macro_node_to_col.insert( std::pair< dof_id_type, unsigned int >(macro_element->node_id(n), 0));
        }
    }
    return;
}

const std::map< dof_id_type, unsigned int >* NodalOverlapUserObject::get_micro_node_to_row() const {
    /*!
    Get the pointer to the micro node to unscaled shape-function matrix column
    */

    return &micro_node_to_row;
}

const std::map< dof_id_type, unsigned int >* NodalOverlapUserObject::get_macro_node_to_col() const {
    /*!
    Get the pointer to the macro node to unscaled shape-function matrix row
    */

    return &macro_node_to_col;
}

void NodalOverlapUserObject::checkElementDepth( const dof_id_type &macro_element, std::map< dof_id_type, int > &element_depths){
    /*!
    Check the depth that the element is in the overlap coupling domain

    Note: This may not be a totally general approach and may give unexpected results for totally submerged macro-scale 
    meshes among others.

    :param dof_id_type macro_element: The id number of the element to check
    :param std::map< dof_id_type, unsigned int >: The depth of the elements explored
    */

    std::map<dof_id_type, std::vector< dof_id_type > >::iterator it;
    std::map<dof_id_type, std::vector< dof_id_type > >::iterator it1;
    std::map< dof_id_type, int >::iterator it2;
    std::vector< int > neighbor_depths;
    int min_pos_depth = 0;
//    std::cout << "macro_element: " << macro_element << "\n";

    it2 = element_depths.find(macro_element); //Check if the macro-element is in the element_depths map
    if (it2 == element_depths.end()){

        element_depths.insert(std::pair< dof_id_type, unsigned int>(macro_element, -1));

        //Get the element's neighbors
        it = macro_element_neighbors.find(macro_element);
        if (it == macro_element_neighbors.end()){
            std::cout << "Error: Element not found in neighbor map\n";
            assert(1==0);
        }

        neighbor_depths.reserve(it->second.size());

        //Iterate through the macro element's neighbors
        for (unsigned int i=0; i<it->second.size(); i++){
//            std::cout << "i: " << i << "\n";
            //Check if the neighbor is not in the overlapped elements (macro_element is on the boundary)
            it1 = macro_element_to_micro_nodes.find(it->second[i]);
            if (it1 == macro_element_to_micro_nodes.end()){
                element_depths[macro_element] = 0; //The macro_element is on the border
                return;
            }
        }

        //Iterate through the macro-element's neighbors now knowing that it isn't on the border.
        for (unsigned int i=0; i<it->second.size(); i++){
            //Because the neighbor is overlapped, check if it has been explored
            it2 = element_depths.find(it->second[i]);
//            std::cout << "neighbor[" << i << "]: " << it->second[i] << "\n";
            if (it2 == element_depths.end()){
//                std::cout << "Exploring neighbor\n";
                checkElementDepth(it2->first, element_depths); //Explore the neighbor
            }

            //Check for any elements which have a depth greater than or equal to zero
            if (element_depths[it2->first] > -0.5){
//                std::cout << "Neighbor depth: " << element_depths[it2->first] << "\n";
                neighbor_depths.push_back(element_depths[it2->first]);
                min_pos_depth = std::min(min_pos_depth, element_depths[it2->first]);
            }
        }

//        std::cout << "neighbor_depths: "; print_vector(neighbor_depths);
//        std::cout << "min_pos_depth: " << min_pos_depth << "\n";
        if (neighbor_depths.size() == 0){ //If all surrounding elements are unknown set the current element to two (not general)
            element_depths[macro_element] = 2;
        }
        else{ //Increment the element depth
            element_depths[macro_element] = min_pos_depth + 1;
        }
//        std::cout << "Element depth: " << element_depths[macro_element] << "\n";
    }
    else if (it2->second < 0){
        std::cout << "Error: Element has a negative depth.\n";
        assert(1==0);
    }
    return;
}

template< class myType >
void
NodalOverlapUserObject::print_vector(std::vector< myType > &v){
    /*!
    print a vector to the terminal

    :param std::vector< myType > v: The vector to be printed.
    */

    for (unsigned int i=0; i<v.size(); i++){
        std::cout << v[i] << " ";
    }
    std::cout << "\n";
}

template< class myType >
void
NodalOverlapUserObject::print_matrix(std::vector< std::vector< myType > > &m){
    /*!
    print a matrix to the terminal

    :param std::vector< std::vector< myType > > m: The matrix to be printed
    */

    for (unsigned int i=0; i<m.size(); i++){
        print_vector(m[i]);
    }
}

void
NodalOverlapUserObject::get_node_info(unsigned int &_num_macro_ghost, unsigned int &_num_macro_free,
                                      unsigned int &_num_micro_ghost, unsigned int &_num_micro_free) const{
    /*!
    Return quantities of interest about the nodes contained in the overlap domain.

    :param unsigned int _num_macro_ghost: The number of ghost nodes in the macro-domain
    :param unsigned int _num_macro_free: The number of free nodes in the macro-domain
    :param unsigned int _num_micro_ghost: The number of ghost nodes in the micro-domain
    :param unsigned int _num_micro_free: The number of free nodes in the micro-domain
    */
    _num_macro_ghost = num_macro_ghost;
    _num_macro_free = num_macro_free;
    _num_micro_ghost = num_micro_ghost;
    _num_micro_free = num_micro_free;
    
}

unsigned int
NodalOverlapUserObject::get_node_elcount(const dof_id_type & id) const{
    /*!
    Return the number of elements the current node is inside. Only greater than 1 for nodes on the boundary of two elements.
    */

    auto it = micro_node_elcount.find(id);
    if (it != micro_node_elcount.end()){
        return it->second;
    }
    else{
        return 1;
    }
}

bool
NodalOverlapUserObject::share_ghost_free_boundary_nodes() const{
    /*!
    Return a boolean indicating if nodes on the boundary between free and ghost macro elements should be shared
    */
    return share_ghost_boundary_nodes;
}

const std::map< dof_id_type, unsigned int>*
NodalOverlapUserObject::get_micro_node_elcount() const{
    /*!
    Return the map from the micro nodes on the boundary to the number of elements they touch.
    */

    return &micro_node_elcount;
}

bool
NodalOverlapUserObject::is_macro_elem_ghost(dof_id_type elnum) const{
    /*!
    Return whether a macro element is a ghost or not.
    */

    auto it = element_depth.find(elnum);
    if (it == element_depth.end()){
        return false;
    }
    else{
        return (it->second >= ghost_depth);
    }
}

//void
//NodalOverlapUserObject::updateMicroElementGPCoords(dof_id_type macro_elem_id, dof_id_type micro_elem_id)
//{
//    /*!
//    Update the positions of the gauss points of a micro element in the local coordinates of the macro element
//    */
//
//    auto _macro_elem = _mesh.elemPtr( macro_elem_id);
//    auto _micro_elem = _mesh.elemPtr( micro_elem_id);
//
//    std::map< dof_id_type, std::vector< Point > > inner;
//    inner.insert(std::make_pair( micro_elem_id, {}));
//
//    //!Compute the local coordinates of the gauss points 
//    std::cout << _assembly << "\n";
//
//    //!Add the macro element to the map if it isn't there
//    std::map< dof_id_type, std::map< dof_id_type, std::vector< Point > > >::iterator it;
//    it = local_micro_gps.find(macro_elem_id);
//    if (it == local_micro_gps.end()){
//        local_micro_gps.insert( std::make_pair(macro_elem_id, inner ));
//    }
//    else {
//        local_micro_gps[macro_elem_id].insert( inner );
//    }
//}
