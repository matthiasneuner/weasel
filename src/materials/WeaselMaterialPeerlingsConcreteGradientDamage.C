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

#include "WeaselMaterialPeerlingsConcreteGradientDamage.h"
#include "WeaselContinuumMechanics.h"
#include <Eigen/Dense>

registerMooseObject("WeaselApp", WeaselMaterialPeerlingsConcreteGradientDamage);

InputParameters
WeaselMaterialPeerlingsConcreteGradientDamage::validParams()
{
  InputParameters params = WeaselMaterialPeerlingsConcreteDamage::validParams();
  params.addClassDescription(
      "Compute stress using a hypoelastic material model from MarmotUserLibrary");
  params.addRequiredParam<Real>("nonlocal_radius", "Nonlocal radius");
  params.addRequiredCoupledVar("nonlocal_damage", "The nonlocal damage variable");
  return params;
}

WeaselMaterialPeerlingsConcreteGradientDamage::WeaselMaterialPeerlingsConcreteGradientDamage(
    const InputParameters & parameters)
  : WeaselMaterialPeerlingsConcreteDamage(parameters),
    _the_nonlocal_radius(getParam<Real>("nonlocal_radius")),
    _equivalent_strain_nonlocal(coupledValue("nonlocal_damage")),
    _equivalent_strain_local(declareProperty<Real>("local_damage")),
    _dequivalent_strain_local_dstrain(
        declarePropertyDerivative<RankTwoTensor>(_base_name + "local_damage", "strain")),
    _dequivalent_strain_local_dequivalent_strain_nonlocal(
        declarePropertyDerivative<Real>(_base_name + "local_damage", "nonlocal_damage")),
    _dstress_dequivalent_strain_nonlocal(
        declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", "nonlocal_damage")),
    _nonlocal_radius(declareProperty<Real>("nonlocal_radius"))
{
}

void
WeaselMaterialPeerlingsConcreteGradientDamage::computeQpProperties()
{
  // Get voigt vectors from MOOSE tensors
  const auto strain = Weasel::voigtFromRankTwoTensor(_strain[_qp], true);
  const auto equivalentStrainNonlocal = _equivalent_strain_nonlocal[_qp];

  // get old kappa
  const auto & kappaOld = _kappa_old[_qp];

  const auto [equivalentStrainLocal, dEquivalentStrainLocal_dStrain] =
      computeModifiedVonMisesStrain(strain);

  const auto [kappaNew, dKappaNew_dEquivalentStrainNonlocal] =
      computeKappa(equivalentStrainNonlocal, kappaOld);

  const auto [omega, dOmega_dKappa] = computeDamage(kappaNew);

  Vector6r stress = (1 - omega) * _Cel * strain;

  Matrix6r dStress_dStrain = (1 - omega) * _Cel;

  Vector6r dStress_dEquivalentStrainNonlocal =
      stress / (1 - omega) * (-dOmega_dKappa * dKappaNew_dEquivalentStrainNonlocal);

  // Convert voigt notation back to MOOSE tensors
  _stress[_qp] = Weasel::rankTwoTensorFromVoigt(stress, false);
  _dstress_dstrain[_qp] = Weasel::rankFourTensorFromVoigt(dStress_dStrain);
  _dstress_dequivalent_strain_nonlocal[_qp] =
      Weasel::rankTwoTensorFromVoigt(dStress_dEquivalentStrainNonlocal, false);

  _equivalent_strain_local[_qp] = equivalentStrainLocal;
  _dequivalent_strain_local_dstrain[_qp] =
      Weasel::rankTwoTensorFromVoigt(dEquivalentStrainLocal_dStrain, false);
  _dequivalent_strain_local_dequivalent_strain_nonlocal[_qp] = 0.0;

  _nonlocal_radius[_qp] = _the_nonlocal_radius;

  // store kappa
  _kappa[_qp] = kappaNew;

  // store damage
  _omega[_qp] = omega;
}
