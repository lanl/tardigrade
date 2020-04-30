//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef DNSFVECTORS_H
#define DNSFVECTORS_H

// MOOSE includes
#include "ElementUserObject.h"
#include "MooseVariableInterface.h"
#include "NodalOverlapUserObject.h"
#include "RankTwoTensor.h"
#include "libmesh/fe.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/bounding_box.h"
#include "libmesh/fe_interface.h"

// Forward Declarations
class DNSFVectors;

template <>
InputParameters validParams<DNSFVectors>();

class DNSFVectors : public ElementUserObject,
                    public MooseVariableInterface<Real>{
    public:
        DNSFVectors(const InputParameters & parameters);

        virtual void initialize() override;
        virtual void execute() override;
        virtual void threadJoin(const UserObject & y) override;
        virtual void finalize() override;


    protected:
        //!The weight
        const MooseArray<Real> & _JxW;
        //!The scaling factor for different coordinate systems
        const MooseArray<Real> & _coord;
        //!The dummy variable object
        MooseVariable & _var;
        //!The test function
        const VariableTestValue & _test;
        const VariableTestGradient & _grad_test;

        //!The nodal overlap user object
        const NodalOverlapUserObject &_nodal_overlap;

        const std::string &_base_name;

        const MaterialProperty<RankTwoTensor> &_stress;

        //!The current quadrature point
        unsigned int _qp;
        //!The test function number
        unsigned int _i;
        //!The dimension being computed
        unsigned int _j;

        //Force Vectors
        std::vector< double > FintQh;

        //Settings
        unsigned int num_micro_dof = 3; //!The number of degrees of freedom for the DNS

        //Overlap information
        unsigned int num_macro_ghost, num_macro_free, num_micro_ghost, num_micro_free;

        //Maps
        const std::map< dof_id_type, unsigned int >* micro_node_to_row;

        Real computeFintQp();
        void computeFVectors();

};

#endif
