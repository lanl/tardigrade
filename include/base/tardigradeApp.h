//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef TARDIGRADEAPP_H
#define TARDIGRADEAPP_H

#include "MooseApp.h"

class tardigradeApp;

class tardigradeApp : public MooseApp
{
public:
  tardigradeApp(InputParameters parameters);
  virtual ~tardigradeApp();
  
  static InputParameters validParams();

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
};

#endif /* TARDIGRADEAPP_H */
