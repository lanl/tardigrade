//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ELEMENTINTEGRATEUSEROBJECT_H
#define ELEMENTINTEGRATEUSEROBJECT_H

// MOOSE includes
#include "ElementUserObject.h"
#include "MooseVariableInterface.h"
//#include "NodalOverlapUserObject.h"
//#include "libmesh/fe.h"
//#include "libmesh/mesh_base.h"
//#include "libmesh/mesh_tools.h"
//#include "libmesh/bounding_box.h"
//#include "libmesh/fe_interface.h"

// Forward Declarations
class ElementIntegrateUserObject;

class ElementIntegrateUserObject : public ElementUserObject,
                                 public MooseVariableInterface<Real>{
    public:
        ElementIntegrateUserObject(const InputParameters & parameters);

        static InputParameters validParams();

        virtual void initialize() override;
        virtual void execute() override;
        virtual void threadJoin(const UserObject & y) override;
        virtual void finalize() override;

        //!Return the nodal volume
        int get_nodal_volume(dof_id_type micro_node_id, double &volume) const;
        double get_nodal_volume(dof_id_type micro_node_id) const;

        //!Compute the nodal density and return it
        int get_nodal_density(dof_id_type micro_node_id, double &density) const;
        double get_nodal_density(dof_id_type micro_node_id) const;

    protected:
        //!The dummy variable object
        MooseVariable & _var;
        //!The test function
        const VariableTestValue & _test;

        //!The current quadrature point
        unsigned int _qp;
        //!The test function number
        unsigned int _i;

        //!The volumes integrated over the surrounding elements
        std::map< dof_id_type, double > integrated_volumes;

        //!The weights integrated over the surrounding elements
        std::map< dof_id_type, double > integrated_weights;

        //!The weighted integral of the densities at the nodes
        std::map< dof_id_type, double > integrated_weighted_densities;

        Real computeQpIntegral();
        void computeIntegral();

        const MaterialProperty<Real> & _density;

};

#endif
