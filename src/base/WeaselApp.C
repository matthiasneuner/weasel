#include "WeaselApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
WeaselApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  params.set<bool>("use_legacy_dirichlet_bc") = false;

  // Do not use legacy material output, i.e., output properties on INITIAL as well as TIMESTEP_END
  params.set<bool>("use_legacy_material_output") = false;

  return params;
}

WeaselApp::WeaselApp(InputParameters parameters) : MooseApp(parameters)
{
  WeaselApp::registerAll(_factory, _action_factory, _syntax);
}

WeaselApp::~WeaselApp() {}

void
WeaselApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"WeaselApp"});
  Registry::registerActionsTo(af, {"WeaselApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
WeaselApp::registerApps()
{
  registerApp(WeaselApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
WeaselApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  WeaselApp::registerAll(f, af, s);
}
extern "C" void
WeaselApp__registerApps()
{
  WeaselApp::registerApps();
}
