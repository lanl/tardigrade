//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef PROJECTIONUSEROBJECT_H
#define PROJECTIONUSEROBJECT_H

// MOOSE includes
#include "NodalUserObject.h"

// Forward Declarations
class NodalOverlapUserObject;

template <>
InputParameters validParams<NodalOverlapUserObject>();

class NodalOverlapUserObject : public NodalUserObject{
public:
  NodalOverlapUserObject(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override {}

protected:

    //!Parameters
    SubdomainName _macroscale_domain; //! The name of the macro-scale subdomain

//    MooseMesh _projection_mesh; //! The mesh to be used in the projection
    SubdomainID _macro_id; //! The id of the macro-scale subdomain
    BoundingBox _macro_bounding_box; //! The bounding box of the macro-scale subdomain
    
};

#endif
