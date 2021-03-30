/* ---------------------------------------------------------------------
 *                              _
 * __      _____  __ _ ___  ___| |
 * \ \ /\ / / _ \/ _` / __|/ _ \ |
 *  \ V  V /  __/ (_| \__ \  __/ |
 *   \_/\_/ \___|\__,_|___/\___|_|
 *
 * Weasel - a simple solid mechanics MOOSE app for teaching purposes
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2021 - today
 *
 * Matthias Neuner matthias.neuner@uibk.ac.at
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of weasel.
 * ---------------------------------------------------------------------
 */

#include "WeaselGradientDamage.h"
#include "DerivativeMaterialInterface.h"

registerMooseObject("WeaselApp", WeaselGradientDamage);

InputParameters
WeaselGradientDamage::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Divergence of stress tensor");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");
  params.addRequiredCoupledVar("nonlocal_damage", "The nonlocal damage field");
  return params;
}

WeaselGradientDamage::WeaselGradientDamage(const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _local_damage(getMaterialPropertyByName<Real>(_base_name + "local_damage")),
    _nonlocal_radius(getMaterialPropertyByName<Real>(_base_name + "nonlocal_radius")),
    _dlocal_damage_dstrain(
        getMaterialPropertyDerivative<RankTwoTensor>(_base_name + "local_damage", "strain")),
    _dlocal_damage_dnonlocal_damage(
        getMaterialPropertyDerivative<Real>(_base_name + "local_damage", "nonlocal_damage")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp)
{
  if (_ndisp != 3)
    mooseError("Weasel kernels are implemented only for 3D!");

  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
WeaselGradientDamage::computeQpResidual()
{
  const auto & NA = _test[_i][_qp];
  const auto & dNA_dx = _grad_test[_i][_qp];
  const auto & dK_dx = _grad_u[_qp];
  const auto & K = _u[_qp];
  const auto & KLocal = _local_damage[_qp];

  const Real lxl = std::pow(_nonlocal_radius[_qp], 2);

  return lxl * (dNA_dx * dK_dx) + NA * (K - KLocal);
}

Real
WeaselGradientDamage::computeQpJacobian()
{
  return computeQpJacobianNonlocalDamage();
}

Real
WeaselGradientDamage::computeQpOffDiagJacobian(unsigned int jvar)
{
  for (unsigned int j = 0; j < _ndisp; ++j)
    if (jvar == _disp_var[j])
      return computeQpJacobianDisplacement(j);

  mooseError("Jacobian for unknown variable requested");
  return 0.0;
}

Real
WeaselGradientDamage::computeQpJacobianDisplacement(unsigned int comp_u)
{
  const auto & NA = _test[_i][_qp];
  const auto & dNB_dx = _grad_phi[_j][_qp];
  const auto & j = comp_u;
  const auto & KLocal = _local_damage[_qp];
  const auto & dKLocal_dE = _dlocal_damage_dstrain[_qp];

  return NA * (0 - 1) * dKLocal_dE.column(j) * dNB_dx;
}

Real
WeaselGradientDamage::computeQpJacobianNonlocalDamage()
{
  const auto & NA = _test[_i][_qp];
  const auto & NB = _phi[_j][_qp];
  const auto & dNA_dx = _grad_test[_i][_qp];
  const auto & dNB_dx = _grad_phi[_j][_qp];

  const Real lxl = std::pow(_nonlocal_radius[_qp], 2);

  return lxl * (dNA_dx * dNB_dx) + NA * NB;
}
