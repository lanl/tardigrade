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
    std::map< dof_id_type, std::vector< dof_id_type > >::iterator it;
    std::vector< dof_id_type >::iterator vec_it;

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
                
//                Point local_point = FEInterface::inverse_map(_macro_element->dim(),
//                                                             libMesh::FEType(_macro_element->default_order()),
//                                                             _macro_element,
//                                                             *_current_point);

                //!Associate the micro-node with a macro-scale element
                it = macro_element_to_micro_nodes.find(_macro_element->id());
                if (it != macro_element_to_micro_nodes.end()){
                    it->second.push_back(_current_node->id());
                }
                else{
                    macro_element_to_micro_nodes.insert( std::pair< dof_id_type, std::vector< dof_id_type > >(_macro_element->id(),
                                                                                                              { _current_node->id()}));
                }

                //!Find the micro-elements associated with a given micro-node
                auto node_elements = micro_node_to_micro_elements.find(_current_node->id());
                for (auto elem_id = node_elements->second.begin(); 
                     elem_id != node_elements->second.end();
                     elem_id++){

                    //!Search for the element id 
                    vec_it = std::find( micro_elements.begin(), micro_elements.end(), *elem_id);
                    
                    //!Store the ids of the elements which are attached to _current_node but haven't been identified yet
                    if (vec_it == micro_elements.end()){
                        micro_elements.push_back(*elem_id);
                    }
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
    for (auto elem_id = pps.micro_elements.begin(); elem_id != pps.micro_elements.end(); elem_id++){

        auto elem_location = find(micro_elements.begin(), micro_elements.end(), *elem_id);
        if (elem_location == micro_elements.end()){
            micro_elements.push_back(*elem_id);
        }
    }
}

void
NodalOverlapUserObject::finalize()
{

    std::cout << "Finalizing Nodal Overlap\n";
    std::cout << "\tNumber of macro-scale elements overlapping microscale domain: " << macro_element_to_micro_nodes.size() << "\n";
    std::cout << "\tNumber of micro-elements associated with overlap:             " << micro_elements.size() << "\n";

}
