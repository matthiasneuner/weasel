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

#pragma once

#include "DerivativeMaterialInterface.h"
#include "WeaselContinuumMechanics.h"

class WeaselMaterialDamagedVonMises : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  WeaselMaterialDamagedVonMises(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// check if the yield condition is violated
  bool checkIfYielding(const Vector6r & trialStress);

  /// simple von Mises yield criterion
  Real yieldFunction(const Vector6r & stress);
  Vector6r dYieldFunction_dStress(const Vector6r & stress);

  /// simple function for computing the increment of equivalent_plastic_strain
  Real deltaKappa(const Vector6r & dStrainPlastic);

  /// assemble the nonlinear equation system rhs of the return mapping algorithm
  VectorXr F(const VectorXr & X);
  /// compute the Jacobian (numerically) using finite differences
  MatrixXr dFdX(const VectorXr & X);

  /// compute scalar isotropic damage based on kappa, softening modulus epsilon_f and an exponential softening law
  std::pair<Real, Real> computeDamage(Real kappa);

  struct ReturnMappingResult
  {
    Vector6r stress;
    Matrix6r dStress_dStrain;

    Real kappa;
    Vector6r dKappa_dStrain;
  };
  /// perform the return mapping algorithm
  ReturnMappingResult performReturnMapping(const Vector6r & trialStress, Real trialKappa);

  /// base name for multiple materials on a block
  const std::string _base_name;

  /// Young's modulus
  const Real _E;
  /// Poisson ratio
  const Real _nu;
  /// Yield strength
  const Real _fcy;
  /// Softening modulus
  const Real _epsilon_f;
  /// elastic stiffness
  const Matrix6r _Cel;

  const MaterialProperty<Real> & _equivalent_plastic_strain_old;
  MaterialProperty<Real> & _equivalent_plastic_strain;

  const MaterialProperty<RankTwoTensor> & _dstrain;
  MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankTwoTensor> & _stress_old;
  MaterialProperty<RankFourTensor> & _dstress_dstrain;
};
