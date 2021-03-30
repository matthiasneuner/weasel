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

#include "WeaselStressDivergence.h"
#include "DerivativeMaterialInterface.h"

registerMooseObject("WeaselApp", WeaselStressDivergence);

InputParameters
WeaselStressDivergence::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Divergence of stress tensor");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction "
                                        "the variable this kernel acts in. (0 for x, "
                                        "1 for y, 2 for z)");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");
  params.addCoupledVar("nonlocal_damage", "The nonlocal damage field");
  return params;
}

WeaselStressDivergence::WeaselStressDivergence(const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _dstress_dstrain(getMaterialPropertyByName<RankFourTensor>(_base_name + "Jacobian_mult")),
    _has_nonlocal_damage(isCoupled("nonlocal_damage")),
    _dstress_dk(
        getMaterialPropertyDerivative<RankTwoTensor>(_base_name + "stress", "nonlocal_damage")),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp),
    _nonlocal_damage_var(_has_nonlocal_damage ? coupled("nonlocal_damage") : 0)
{
  if (_ndisp != 3)
    mooseError("Weasel kernels are implemented only for 3D!");

  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
WeaselStressDivergence::computeQpResidual()
{
  // We abandon MOOSE's 'creative' naming scheme here:
  //
  // Node index A = _i
  // direction test function i = _component
  //
  // What we need to compute
  //
  //               ∂ N_A
  // R_(Ai)     =  ----- * σₖᵢ
  //               ∂ xₖ

  const auto & dNA_dx = _grad_test[_i][_qp];
  const auto & i = _component;
  const auto & S = _stress[_qp];

  return dNA_dx * S.column(i);
}

Real
WeaselStressDivergence::computeQpJacobian()
{
  const auto ivar = _var.number();

  for (unsigned int i = 0; i < _ndisp; ++i)
    if (ivar == _disp_var[i])
      return computeQpJacobianDisplacement(_component, _component);

  mooseError("Jacobian for unknown variable requested");
  return 0.0;
}

Real
WeaselStressDivergence::computeQpOffDiagJacobian(unsigned int jvar)
{
  for (unsigned int j = 0; j < _ndisp; ++j)
    if (jvar == _disp_var[j])
      return computeQpJacobianDisplacement(_component, j);

  if (_has_nonlocal_damage && jvar == _nonlocal_damage_var)
    return computeQpJacobianNonlocalDamage(_component);

  mooseError("Jacobian for unknown variable requested");
  return 0.0;
}

Real
WeaselStressDivergence::computeQpJacobianDisplacement(unsigned int comp_i, unsigned int comp_j)
{
  // We abandon MOOSE's 'creative' naming scheme here:
  //
  // Node index A = _i
  // Node index B = _j
  // direction test function i = comp_i
  // direction trial function j = comp_j
  //
  // What we need to compute
  //
  //              ∂ f_Aᵢ     ∂ N_A   ∂ σₖᵢ     ∂ Δε ₗₘ
  // J_(Ai)(Bj) = -------- = ----- * ------- * -------
  //              ∂ Δqᵘ_Bⱼ    ∂ xₖ   ∂ Δε ₗₘ   ∂ Δqᵘ_Bⱼ
  //
  //
  //              ∂ f_Aᵢ     ∂ N_A   ∂ σₖᵢ     ∂ N_B
  //            = -------- = ----- * ------  * -----
  //              ∂ Δqᵘ_Bⱼ    ∂ xₖ   ∂ Δε ₗⱼ   ∂ xₗ

  const auto & dNA_dx = _grad_test[_i][_qp];
  const auto & dNB_dx = _grad_phi[_j][_qp];
  const auto & i = comp_i;
  const auto & j = comp_j;
  const auto & dS_dE = _dstress_dstrain[_qp];

  // Moose tensors offers no higher horder contraction possibilites, so we need to do it by hand
  Real dResidual_Ai_dNodeDisplacement_Bj = 0;
  for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++)
      dResidual_Ai_dNodeDisplacement_Bj += dNA_dx(k) * dS_dE(k, i, l, j) * dNB_dx(l);

  return dResidual_Ai_dNodeDisplacement_Bj;
}

Real
WeaselStressDivergence::computeQpJacobianNonlocalDamage(unsigned int comp_i)
{
  // We abandon MOOSE's 'creative' naming scheme here:
  //
  // Node index A = _i
  // Node index B = _j
  // direction test function i = comp_i
  //
  // What we need to compute
  //
  //              ∂ f_Aᵢ     ∂ N_A   ∂ σₖᵢ
  // J_(Ai)(Bj) = -------- = ----- * ------- * N_B
  //              ∂ Δqᵘ_Bⱼ    ∂ xₖ   ∂ K
  //

  const auto & dNA_dx = _grad_test[_i][_qp];
  const auto & i = comp_i;
  const auto & NB = _phi[_j][_qp];
  const auto & dS_dK = _dstress_dk[_qp];

  return dNA_dx * dS_dK.column(i) * NB;
}
