//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "tardigradeTestApp.h"
#include "tardigradeApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
tardigradeTestApp::validParams()
{
  InputParameters params = tardigradeApp::validParams();
  return params;
}

tardigradeTestApp::tardigradeTestApp(InputParameters parameters) : MooseApp(parameters)
{
  tardigradeTestApp::registerAll(
     _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

tardigradeTestApp::~tardigradeTestApp() {}

void
tardigradeTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  tardigradeApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"tardigradeTestApp"});
    Registry::registerActionsTo(af, {"tardigradeTestApp"});
  }
}

void
tardigradeTestApp::registerApps()
{
  registerApp(tardigradeApp);
  registerApp(tardigradeTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
tardigradeTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  tardigradeTestApp::registerAll(f, af, s);
}
extern "C" void
tardigradeTestApp__registerApps()
{
  tardigradeTestApp::registerApps();
}
