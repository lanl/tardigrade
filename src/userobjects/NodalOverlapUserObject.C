//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NodalOverlapUserObject.h"

#include "libmesh/bounding_box.h"
#include "libmesh/mesh_tools.h"

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

                //!Store the node in the id to shape-function col index map
                updateMicroNodeColMap(_current_node->id());                

                //!Store the macro-element's nodes in the id to shape-function row map
                updateMacroNodeRowMap(_macro_element);

                //!Find the micro-elements associated with a given micro-node
                auto node_elements = micro_node_to_micro_elements.find(_current_node->id());
                for (auto elem_id = node_elements->second.begin(); 
                     elem_id != node_elements->second.end();
                     elem_id++){

                    //!Update micro_elements if required
                    updateMicroElements(_macro_element->id(), *elem_id);

                }
                break; //A node can only appear in a single element
            }
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
    for (it1=macro_element_to_micro_nodes.begin(); it1 != macro_element_to_micro_nodes.end(); it1++){
        std::cout << "\t  Macro element id: " << it1->first << " Number of nodes: " << it1->second.size() << "\n";
    }

    std::cout << "\tNumber of micro-elements associated with overlap:             " << micro_elements.size() << "\n";

    //!Update the micro node to column map
    std::map< dof_id_type, unsigned int >::iterator it;
    unsigned int index=0;

    for (it=micro_node_to_col.begin(); it!=micro_node_to_col.end(); it++){
        it->second = index;
        index++;
    }

    //!Update the macro node to row map
    index = 0;
    for (it=macro_node_to_row.begin(); it!=macro_node_to_row.end(); it++){
        it->second = index;
        index++;
    }

//    std::map< dof_id_type, std::vector< double > >::iterator it;
//    std::cout << "Node local coordinates\n";
//    for (it=NodeLocalCoordinates.begin(); it != NodeLocalCoordinates.end(); it++){
//        std::cout << it->first << ": " << it->second[0] << " " << it->second[1] << " " << it->second[2] << "\n";
//    }

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

void NodalOverlapUserObject::updateMicroNodeColMap(const dof_id_type micro_node_id){
    /*!
    Add the indicated node id to the id to shape-function matrix column map. The column number will be 
    determined at finalize. Note that the number is the order. There will be a scale factor applied due to 
    the number of degrees of freedom to be mapped from the micro to macro scales.


    :param dof_id_type micro_node_id: The node id to add
    */

//    std::map< dof_id_type, unsigned int >::iterator it;
//    it = micro_node_to_col.find(micro_node_id);
//    if (it == micro_node_to_col.end()){
//        micro_node_to_col.insert( std::pair< dof_id_type, unsigned int>(micro_node_id, 0));
//    }

    micro_node_to_col.insert( std::pair< dof_id_type, unsigned int>(micro_node_id, 0));
    return;
}

void NodalOverlapUserObject::updateMacroNodeRowMap(const Elem* macro_element){
    /*!
    Add the nodes of the macro element to the id to shape-function matrix column map. The row number will 
    be determined at finalize. Note that the number is the order. There will be a scale factor applied due to 
    the number of degrees of freedom to be mapped from the micro to macro scales.

    :param Elem* macro_element: The pointer to the macro-scale element.
    */

    std::map< dof_id_type, unsigned int >::iterator it;

    for (unsigned int n=0; n<macro_element->n_nodes(); n++){
        it = macro_node_to_row.find(macro_element->node_id(n));
        if (it == macro_node_to_row.end()){
            macro_node_to_row.insert( std::pair< dof_id_type, unsigned int >(macro_element->node_id(n), 0));
        }
    }
    return;
}

const std::map< dof_id_type, unsigned int >* NodalOverlapUserObject::get_micro_node_to_col() const {
    /*!
    Get the pointer to the micro node to unscaled shape-function matrix column
    */

    return &micro_node_to_col;
}

const std::map< dof_id_type, unsigned int >* NodalOverlapUserObject::get_macro_node_to_row() const {
    /*!
    Get the pointer to the macro node to unscaled shape-function matrix row
    */

    return &macro_node_to_row;
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
