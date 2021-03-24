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

#include "WeaselMaterialDamagedVonMises.h"
#include "WeaselContinuumMechanics.h"
#include <Eigen/Dense>

registerMooseObject("WeaselApp", WeaselMaterialDamagedVonMises);

InputParameters
WeaselMaterialDamagedVonMises::validParams()
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
  params.addRequiredParam<Real>("fcy", "Yield stress");
  params.addRequiredParam<Real>("epsilon_f", "Softening modulus");
  return params;
}

WeaselMaterialDamagedVonMises::WeaselMaterialDamagedVonMises(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _E(getParam<Real>("E")),
    _nu(getParam<Real>("nu")),
    _fcy(getParam<Real>("fcy")),
    _epsilon_f(getParam<Real>("epsilon_f")),
    _Cel(Weasel::Elasticity::stiffnessTensor(_E, _nu)),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _equivalent_plastic_strain_old(getMaterialPropertyOld<Real>("equivalent_plastic_strain")),
    _equivalent_plastic_strain(declareProperty<Real>("equivalent_plastic_strain")),
    _dstrain(getMaterialProperty<RankTwoTensor>("strain_increment")),
    _stress(declareProperty<RankTwoTensor>(_base_name + "stress")),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "stress")),
    _dstress_dstrain(declareProperty<RankFourTensor>(_base_name + "Jacobian_mult"))
{
}

void
WeaselMaterialDamagedVonMises::initQpStatefulProperties()
{
  // Initalize stress to zero
  _stress[_qp] = Weasel::rankTwoTensorFromVoigt(Vector6r::Zero(), false);
  // Initialize inelastic strain to zero
  _equivalent_plastic_strain[_qp] = 0.0;
}

void
WeaselMaterialDamagedVonMises::computeQpProperties()
{
  // Get voigt vectors from MOOSE tensors
  const auto strainIncrement = Weasel::voigtFromRankTwoTensor(_dstrain[_qp], true);
  const auto stressOld = Weasel::voigtFromRankTwoTensor(_stress_old[_qp], false);
  Real kappaOld = _equivalent_plastic_strain_old[_qp];

  const auto [omegaOld, dOmegaOld_dKappa] = computeDamage(kappaOld);

  Vector6r stress;
  Matrix6r dStress_dStrain;
  Real kappa;
  Real omega;

  // Compute trial stress
  // σᵗʳₙ₊₁ = σₙ + Cᵉˡ : Δ ε
  const auto trialStress = stressOld + _Cel * strainIncrement;

  // Compute effective trial stress
  const auto effectiveTrialStress = trialStress * 1. / (1 - omegaOld);
  const auto trialKappa = kappaOld;

  if (checkIfYielding(effectiveTrialStress))
  {
    // yielding, perform return mapping
    const auto [effectiveStressNew, dEffectiveStressNew_dStrain, kappaNew, dKappa_dStrain] =
        performReturnMapping(effectiveTrialStress, trialKappa);

    // update of damage
    auto [omega, dOmega_dKappa] = computeDamage(kappaNew);

    // new stress is the back projected stess, tangent stiffness is the damaged elastoplastic
    // stiffness
    stress = effectiveStressNew * (1 - omega);
    dStress_dStrain = dEffectiveStressNew_dStrain * (1 - omega) + effectiveStressNew * (1 - omega) *
                                                                      -1 * dOmega_dKappa *
                                                                      dKappa_dStrain.transpose();

    kappa = kappaNew;
  }
  else
  {
    // no yielding, resulting stress is trial stress, tangent stiffness is elastic stiffness
    stress = trialStress * (1 - omega);
    dStress_dStrain = _Cel * (1 - omega);

    // kappa did not change
    kappa = trialKappa;
  }

  // Convert voigt notation back to MOOSE tensors
  _stress[_qp] = Weasel::rankTwoTensorFromVoigt(stress, false);
  _dstress_dstrain[_qp] = Weasel::rankFourTensorFromVoigt(dStress_dStrain);

  // store inelastic equivalent strain
  _equivalent_plastic_strain[_qp] = kappa;
}

Real
WeaselMaterialDamagedVonMises::yieldFunction(const Vector6r & stress)
{
  const Real J2 = Weasel::Invariants::J2(stress);

  // f(σ) = √{3 J₂} − fₜ
  return std::sqrt(3 * J2) - _fcy;
}

Vector6r
WeaselMaterialDamagedVonMises::dYieldFunction_dStress(const Vector6r & stress)
{

  const Real J2 = Weasel::Invariants::J2(stress);
  const Vector6r dJ2_dStress = Weasel::Invariants::dJ2_dStress(stress);

  // ∂f/∂σ  = −0.5 / √{3 J₂} * 3 * ∂J₂/∂ σ
  return -0.5 / (std::sqrt(3 * J2)) * 3 * dJ2_dStress;
}

Real
WeaselMaterialDamagedVonMises::deltaKappa(const Vector6r & dStrainPlastic)
{
  // Correct shear terms for J2 ( of strain ) computation
  Vector6r dStrainPlasticCorrected;
  dStrainPlasticCorrected << dStrainPlastic(0), dStrainPlastic(1), dStrainPlastic(2),
      dStrainPlastic(3) / 2, dStrainPlastic(4) / 2, dStrainPlastic(5) / 2;

  const Real J2 = Weasel::Invariants::J2(dStrainPlasticCorrected);

  return std::sqrt(3 * J2);
}

bool
WeaselMaterialDamagedVonMises::checkIfYielding(const Vector6r & stress)
{
  const auto f = yieldFunction(stress);
  return f > 1e-12;
}

VectorXr
WeaselMaterialDamagedVonMises::F(const VectorXr & X)
{
  const Vector6r & stress = X.head(6);
  const Real & kappa = X(6);
  const Real & dLambda = X(7);

  VectorXr F(X);

  /* F ( X) =
   * / σₙ₊₁ + Cᵉˡ : Δ ε \
   * | ...              |
   * | ...              |
   * | ...              |
   * | ...              |
   * | ...              |
   * | κₙ₊₁ - Δκ        |
   * \ f(σₙ₊₁)          /
   */
  const Vector6r dStrainPlastic = dLambda * dYieldFunction_dStress(stress);

  // clang-format off
  F.head(6) = stress + _Cel * dStrainPlastic;
  F(6)      = kappa - deltaKappa ( dStrainPlastic ) ; 
  F(7)      = yieldFunction(stress);
  // clang-format on  

  return F;
}

MatrixXr
WeaselMaterialDamagedVonMises::dFdX(const VectorXr & X)
{
  const auto xSize = X.rows();
  MatrixXr J(xSize, xSize);

  VectorXr leftX(xSize);
  VectorXr rightX(xSize);


  // compute Jacobian of the return mapping system using central differences !
  
  for (auto i = 0; i < xSize; i++)
  {
    Real volatile h = std::max(1.0, std::abs(X(i))) * 1e-5;
    // clang-format off
        leftX  = X;
        rightX = X;
        leftX( i ) -= h;
        rightX( i ) += h;

        J.col( i ) = ( ( F( rightX) ) - F( leftX ) ) 
                   / //-----------------------------
                                ( 2 * h );
    // clang-format on
  }

  return J;
}

WeaselMaterialDamagedVonMises::ReturnMappingResult
WeaselMaterialDamagedVonMises::performReturnMapping(const Vector6r & trialStress, Real trialKappa)
{
  /* Solve the equation system in the unknowns σₙ₊₁ and Δλ
   *
   * R =    LHS      -      F ( X )
   * / \ = / σᵗʳₙ₊₁ \   / σₙ₊₁ + Cᵉˡ : Δεᵖ \
   * | | = | ...    |   | ...              |
   * | | = | ...    |   | ...              |
   * | | = | ...    |   | ...              |
   * | | = | ...    |   | ...              |
   * | | = | ...    |   | ...              |
   * | | = | κₙ     |   | κₙ₊₁ - Δκ        |
   * \ / = \ 0      /   \ f(σₙ₊₁)          /
   *
   *                 ∂f(σₙ₊₁)
   * with Δεᵖ = Δλ * --------
   *                 ∂σₙ₊₁
   */

  constexpr auto sizeEq = 6 + 1 + 1;

  VectorXr residual(sizeEq);

  VectorXr LHS(sizeEq);
  // left hand side vector X Layout : LHS  = [ σ₁₁, σ₂₂, σ₃₃, σ₁₂, σ₁₃, σ₂₃, κₙ, 0 ]ᵀ
  LHS.head(6) = trialStress;
  LHS(6) = trialKappa;
  LHS(7) = 0.0;

  // X contains the 7 unknowns
  VectorXr X(sizeEq);
  VectorXr dX = VectorXr::Zero(sizeEq);

  // solution vector X Layout : X = [ σ₁₁, σ₂₂, σ₃₃, σ₁₂, σ₁₃, σ₂₃, κₙ₊₁, Δλ ]
  // intial guess for the unknowns: assume the trial state!
  X = LHS;
  residual = LHS - F(X);
  auto resNorm = residual.norm();

  int iterationCounter = 0;
  while (resNorm > 1e-12 && iterationCounter < 10)
  {

    // ΔX = J⁻¹ * R   =  = (∂F(Xₙ₊₁)/∂ Xₙ₊₁)⁻¹ * R
    dX = dFdX(X).colPivHouseholderQr().solve(residual);

    X += dX;
    residual = LHS - F(X);
    resNorm = residual.norm();

    iterationCounter++;
  }

  if (resNorm > 1e-12)
    throw MooseException("Weasel von Mises material failed in return mapping!");

  // new converged stress are the first 6 entries of solution vector X
  const Vector6r newStress = X.head(6);

  // new converged kappa is the 7th entry
  const Real newKappa = X(6);

  /* Compute elasto plastic stiffness according to following scheme
   *
   * ∂ σₙ₊₁   ∂σₙ₊₁      ∂σᵗʳₙ₊₁
   * ------ = -------  *  -------
   * ∂ Δ ε    ∂σᵗʳₙ₊₁    ∂ Δ ε
   *
   * with
   *          ∂σᵗʳₙ₊₁
   *          -------  = Cᵉˡ
   *          ∂ Δ ε
   *
   * and where the relation F( X ) = LHS for the converged sate is exploited:
   *
   * ∂σₙ₊₁               / ∂ F(Xₙ₊₁) \
   * -------  = inverse  | --------- |
   * ∂σᵗʳₙ₊₁             \ ∂ Xₙ₊₁    /
   *                                  [ 6x6 topleft ]
   *
   */

  Eigen::Matrix<Real, sizeEq, 6> dLHS_dStrain;

  dLHS_dStrain.block<6, 6>(0, 0) = _Cel;
  dLHS_dStrain.row(6) = Vector6r::Zero();
  dLHS_dStrain.row(7) = Vector6r::Zero();

  MatrixXr dX_dStrain = dFdX(X).colPivHouseholderQr().solve(dLHS_dStrain);

  const Matrix6r dNewStress_dStrain = dX_dStrain.block<6, 6>(0, 0);
  const Vector6r dNewKappa_dStrain = dX_dStrain.row(6);

  return {.stress = newStress,
          .dStress_dStrain = dNewStress_dStrain,

          .kappa = newKappa,
          .dKappa_dStrain = dNewKappa_dStrain};
}

std::pair<Real, Real>
WeaselMaterialDamagedVonMises::computeDamage(Real kappa)
{
  const Real omega = 1 - std::exp(-kappa / _epsilon_f);
  const Real dOmega_dKappa = 0 - std::exp(-kappa / _epsilon_f) * -1. / _epsilon_f;

  return {omega, dOmega_dKappa};
}
