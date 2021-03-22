//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "WeaselTestApp.h"
#include "WeaselApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
WeaselTestApp::validParams()
{
  InputParameters params = WeaselApp::validParams();
  return params;
}

WeaselTestApp::WeaselTestApp(InputParameters parameters) : MooseApp(parameters)
{
  WeaselTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

WeaselTestApp::~WeaselTestApp() {}

void
WeaselTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  WeaselApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"WeaselTestApp"});
    Registry::registerActionsTo(af, {"WeaselTestApp"});
  }
}

void
WeaselTestApp::registerApps()
{
  registerApp(WeaselApp);
  registerApp(WeaselTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
WeaselTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  WeaselTestApp::registerAll(f, af, s);
}
extern "C" void
WeaselTestApp__registerApps()
{
  WeaselTestApp::registerApps();
}
