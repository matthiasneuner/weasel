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

#include "WeaselMaterialLinearElastic.h"
#include "WeaselContinuumMechanics.h"

registerMooseObject("WeaselApp", WeaselMaterialLinearElastic);

InputParameters
WeaselMaterialLinearElastic::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Compute stress using a hypoelastic material model from MarmotUserLibrary");
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  params.addRequiredParam<Real>("E", "Young's modulus");
  params.addRequiredParam<Real>("nu", "Poisson ratio");
  return params;
}

WeaselMaterialLinearElastic::WeaselMaterialLinearElastic(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _E(getParam<Real>("E")),
    _nu(getParam<Real>("nu")),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _dstrain(getMaterialProperty<RankTwoTensor>("strain_increment")),
    _stress(declareProperty<RankTwoTensor>(_base_name + "stress")),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "stress")),
    _dstress_dstrain(declareProperty<RankFourTensor>(_base_name + "Jacobian_mult"))
{
}

void
WeaselMaterialLinearElastic::initQpStatefulProperties()
{
  // Initialize stress to zero
  _stress[_qp] = Weasel::rankTwoTensorFromVoigt(Vector6r::Zero(), false);
}

void
WeaselMaterialLinearElastic::computeQpProperties()
{
  // Get voigt vectors from MOOSE tensors
  const auto strainIncrement = Weasel::voigtFromRankTwoTensor(_dstrain[_qp], true);
  const auto stressOld = Weasel::voigtFromRankTwoTensor(_stress_old[_qp], false);

  // Compute elastic stiffnes Cᵉˡ from E and ν
  const auto Cel = Weasel::Elasticity::stiffnessTensor(_E, _nu);

  // σₙ₊₁ = σₙ + Cᵉˡ : Δ ε
  const auto stress = stressOld + Cel * strainIncrement;

  // Convert voigt notation back to MOOSE tensors
  _stress[_qp] = Weasel::rankTwoTensorFromVoigt(stress, false);
  _dstress_dstrain[_qp] = Weasel::rankFourTensorFromVoigt(Cel);
}
