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
#include "ElementUserObject.h"

// Forward Declarations
class ProjectionUserObject;

template <>
InputParameters validParams<ProjectionUserObject>();

class ProjectionUserObject : public ElementUserObject{
public:
  ProjectionUserObject(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override {}

protected:
    SubdomainName _macroscale_domain;
};

#endif
