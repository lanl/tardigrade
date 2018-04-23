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

template <>
InputParameters
validParams<tardigradeTestApp>()
{
  InputParameters params = validParams<tardigradeApp>();
  return params;
}

tardigradeTestApp::tardigradeTestApp(InputParameters parameters) : MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  tardigradeApp::registerObjectDepends(_factory);
  tardigradeApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  tardigradeApp::associateSyntaxDepends(_syntax, _action_factory);
  tardigradeApp::associateSyntax(_syntax, _action_factory);

  Moose::registerExecFlags(_factory);
  ModulesApp::registerExecFlags(_factory);
  tardigradeApp::registerExecFlags(_factory);

  bool use_test_objs = getParam<bool>("allow_test_objects");
  if (use_test_objs)
  {
    tardigradeTestApp::registerObjects(_factory);
    tardigradeTestApp::associateSyntax(_syntax, _action_factory);
    tardigradeTestApp::registerExecFlags(_factory);
  }
}

tardigradeTestApp::~tardigradeTestApp() {}

void
tardigradeTestApp::registerApps()
{
  registerApp(tardigradeApp);
  registerApp(tardigradeTestApp);
}

void
tardigradeTestApp::registerObjects(Factory & /*factory*/)
{
  /* Uncomment Factory parameter and register your new test objects here! */
}

void
tardigradeTestApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
  /* Uncomment Syntax and ActionFactory parameters and register your new test objects here! */
}

void
tardigradeTestApp::registerExecFlags(Factory & /*factory*/)
{
  /* Uncomment Factory parameter and register your new execute flags here! */
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
tardigradeTestApp__registerApps()
{
  tardigradeTestApp::registerApps();
}

// External entry point for dynamic object registration
extern "C" void
tardigradeTestApp__registerObjects(Factory & factory)
{
  tardigradeTestApp::registerObjects(factory);
}

// External entry point for dynamic syntax association
extern "C" void
tardigradeTestApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  tardigradeTestApp::associateSyntax(syntax, action_factory);
}

// External entry point for dynamic execute flag loading
extern "C" void
tardigradeTestApp__registerExecFlags(Factory & factory)
{
  tardigradeTestApp::registerExecFlags(factory);
}
