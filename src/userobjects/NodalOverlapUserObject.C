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
                std::vector< double > local_coords(LIBMESH_DIM, 0);
//                for (unsigned int i=0; i<local_coords.size(); i++){local_coords[i] = local_point(i);}
//                NodeLocalCoordinates.insert( std::pair< dof_id_type, Point >(_current_node->id(), local_point));

                //!Associate the micro-node with a macro-scale element
                updateMacroMicroMap(_macro_element->id(), _current_node->id(), local_point);

                //!Find the micro-elements associated with a given micro-node
                auto node_elements = micro_node_to_micro_elements.find(_current_node->id());
                for (auto elem_id = node_elements->second.begin(); 
                     elem_id != node_elements->second.end();
                     elem_id++){

                    //!Update micro_elements if required
                    updateMicroElements(_macro_element->id(), *elem_id);

                }
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
    for (it1=macro_element_to_micro_nodes.begin(); it1 != macro_element_to_micro_nodes.end(); it1++){
        std::cout << "\t  Macro element id: " << it1->first << " Number of nodes: " << it1->second.size() << "\n";
    }

    std::cout << "\tNumber of micro-elements associated with overlap:             " << micro_elements.size() << "\n";
    std::map< dof_id_type, std::vector< double > >::iterator it;
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
    auto it = macro_element_to_micro_positions.find(macro_element_id);
    if (it == macro_element_to_micro_positions.end()){
        return NULL;
    }
    else{
        return &it->second;
    }
    return NULL;
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
