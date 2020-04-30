#include "tardigradeApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

template <>
InputParameters
validParams<tardigradeApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

tardigradeApp::tardigradeApp(InputParameters parameters) : MooseApp(parameters)
{
  tardigradeApp::registerAll(_factory, _action_factory, _syntax);
}

tardigradeApp::~tardigradeApp() {}

void
tardigradeApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"tardigradeApp"});
  Registry::registerActionsTo(af, {"tardigradeApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
tardigradeApp::registerApps()
{
  registerApp(tardigradeApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
tardigradeApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  tardigradeApp::registerAll(f, af, s);
}
extern "C" void
tardigradeApp__registerApps()
{
  tardigradeApp::registerApps();
}
