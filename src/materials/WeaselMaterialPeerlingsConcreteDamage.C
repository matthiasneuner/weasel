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

#include "WeaselMaterialPeerlingsConcreteDamage.h"
#include "WeaselContinuumMechanics.h"
#include <Eigen/Dense>

registerMooseObject("WeaselApp", WeaselMaterialPeerlingsConcreteDamage);

InputParameters
WeaselMaterialPeerlingsConcreteDamage::validParams()
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
  params.addRequiredParam<Real>("kappa_0", "material strength");
  params.addRequiredParam<Real>("alpha", "damage material parameter");
  params.addRequiredParam<Real>("beta", "damage material parameter");
  params.addRequiredParam<Real>("k", "damage material parameter");
  return params;
}

WeaselMaterialPeerlingsConcreteDamage::WeaselMaterialPeerlingsConcreteDamage(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _E(getParam<Real>("E")),
    _nu(getParam<Real>("nu")),
    _kappa_0(getParam<Real>("kappa_0")),
    _alpha(getParam<Real>("alpha")),
    _beta(getParam<Real>("beta")),
    _k(getParam<Real>("k")),
    _Cel(Weasel::Elasticity::stiffnessTensor(_E, _nu)),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _kappa_old(getMaterialPropertyOld<Real>("kappa")),
    _kappa(declareProperty<Real>("kappa")),
    _omega(declareProperty<Real>("damage")),
    _strain(getMaterialProperty<RankTwoTensor>("total_strain")),
    _stress(declareProperty<RankTwoTensor>(_base_name + "stress")),
    _dstress_dstrain(declareProperty<RankFourTensor>(_base_name + "Jacobian_mult"))
{
}

void
WeaselMaterialPeerlingsConcreteDamage::initQpStatefulProperties()
{
  // Initialize stateful kappa to zero
  _kappa[_qp] = 0.0;
}

void
WeaselMaterialPeerlingsConcreteDamage::computeQpProperties()
{

  // Get voigt vectors from MOOSE tensors
  const auto strain = Weasel::voigtFromRankTwoTensor(_strain[_qp], true);

  // get old kappa
  const auto & kappaOld = _kappa_old[_qp];

  // update in equivalent strain
  const auto [equivalentStrain, dEquivalentStrain_dStrain] = computeModifiedVonMisesStrain(strain);

  // update kappa
  const auto [kappaNew, dKappaNew_dEquivalentStrain] = computeKappa(equivalentStrain, kappaOld);

  // compute damage based on new kappa
  const auto [omega, dOmega_dKappa] = computeDamage(kappaNew);

  // update stress and compute algorithmic tanget
  Vector6r stress = (1 - omega) * _Cel * strain;

  // clang-format off
  Matrix6r dStress_dStrain =  stress / ( 1 - omega) * ( - dOmega_dKappa * 
                                                          dKappaNew_dEquivalentStrain * 
                                                          dEquivalentStrain_dStrain.transpose() ) 
                               + ( 1 - omega) * _Cel;
  // clang-format on

  // Convert voigt notation back to MOOSE tensors
  _stress[_qp] = Weasel::rankTwoTensorFromVoigt(stress, false);
  _dstress_dstrain[_qp] = Weasel::rankFourTensorFromVoigt(dStress_dStrain);

  // store kappa
  _kappa[_qp] = kappaNew;

  // make damage available as a material property
  _omega[_qp] = omega;
}

std::pair<Real, Vector6r>
WeaselMaterialPeerlingsConcreteDamage::computeModifiedVonMisesStrain(const Vector6r & strain)
{
  const auto & k = _k;
  const auto & v = _nu;
  const Real I1 = Weasel::Invariants::I1(strain);
  const static Vector6r dI1_dStrain = (Vector6r() << 1, 1, 1, 0, 0, 0).finished();

  const Real J2 = Weasel::Invariants::J2Strain(strain);
  Vector6r dJ2_dStrain = Weasel::Invariants::dJ2Strain_dStrain(strain);

  const Real a = std::pow((k - 1) / (1 - 2 * v), 2);
  const Real b = (k - 1) / (2 * k * (1 - 2 * v));
  const Real c = (2 * k) / ((1 + v) * (1 + v));

  // clang-format off
  const Real aux = std::sqrt ( a *I1*I1  + c * J2  );

  Real dAux_dI1 = 0.0;
  Real dAux_dJ2 = 0.0;

  if ( aux >= 1e-15){
    dAux_dI1 = 0.5 / aux *  a *I1*2      ;
    dAux_dJ2 = 0.5 / aux *              c;
  }
  
  const Real  eq     = b * I1 + 1./(2*k) * aux;
  const Real dEq_dI1 = b      + 1./(2*k) * dAux_dI1;
  const Real dEq_dJ2 =        + 1./(2*k) * dAux_dJ2;
  // clang-format on

  return {eq, dEq_dI1 * dI1_dStrain + dEq_dJ2 * dJ2_dStrain};
}

std::pair<Real, Real>
WeaselMaterialPeerlingsConcreteDamage::computeKappa(Real modifiedVonMisesStrain, Real kappaOld)
{
  const Real kappa = std::max(kappaOld, modifiedVonMisesStrain);

  if (modifiedVonMisesStrain >= kappaOld)
    return {kappa, 1};
  else
    return {kappa, 0};
}

std::pair<Real, Real>
WeaselMaterialPeerlingsConcreteDamage::computeDamage(Real kappa)
{
  if (kappa < _kappa_0)
    return {0, 0};

  // clang-format off
  const Real omega =         1 - (    _kappa_0 / kappa )         *  ( (1 - _alpha)  + _alpha * std::exp( - _beta * ( kappa - _kappa_0 ) ) );

  const Real dOmega_dKappa = 0 - (
                                 ( - _kappa_0 / (kappa*kappa) )  *  ( (1 - _alpha)  + _alpha * std::exp( - _beta * ( kappa - _kappa_0 ) ) ) + 
                                 (   _kappa_0 / kappa )          *  (               + _alpha * std::exp( - _beta * ( kappa - _kappa_0 ) )   *  - _beta ));
  // clang-format on
  return {omega, dOmega_dKappa};
}
