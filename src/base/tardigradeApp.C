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
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  tardigradeApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  tardigradeApp::associateSyntax(_syntax, _action_factory);

  Moose::registerExecFlags(_factory);
  ModulesApp::registerExecFlags(_factory);
  tardigradeApp::registerExecFlags(_factory);
}

tardigradeApp::~tardigradeApp() {}

void
tardigradeApp::registerApps()
{
  registerApp(tardigradeApp);
}

void
tardigradeApp::registerObjects(Factory & factory)
{
    Registry::registerObjectsTo(factory, {"tardigradeApp"});
}

void
tardigradeApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & action_factory)
{
  Registry::registerActionsTo(action_factory, {"tardigradeApp"});

  /* Uncomment Syntax parameter and register your new production objects here! */
}

void
tardigradeApp::registerObjectDepends(Factory & /*factory*/)
{
}

void
tardigradeApp::associateSyntaxDepends(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}

void
tardigradeApp::registerExecFlags(Factory & /*factory*/)
{
  /* Uncomment Factory parameter and register your new execution flags here! */
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
tardigradeApp__registerApps()
{
  tardigradeApp::registerApps();
}

extern "C" void
tardigradeApp__registerObjects(Factory & factory)
{
  tardigradeApp::registerObjects(factory);
}

extern "C" void
tardigradeApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  tardigradeApp::associateSyntax(syntax, action_factory);
}

extern "C" void
tardigradeApp__registerExecFlags(Factory & factory)
{
  tardigradeApp::registerExecFlags(factory);
}
