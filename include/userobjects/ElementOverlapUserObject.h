//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ELEMENTOVERLAPUSEROBJECT_H
#define ELEMENTOVERLAPUSEROBJECT_H

// MOOSE includes
#include "ElementUserObject.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/bounding_box.h"
#include "libmesh/fe_interface.h"

// Forward Declarations
class ElementOverlapUserObject;

template <>
InputParameters validParams<ElementOverlapUserObject>();

class ElementOverlapUserObject : public ElementUserObject{
public:
  ElementOverlapUserObject(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override;

};

#endif
