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
/**
 * Implementation of the (local) damage mode for concrete after
 * Peerlings, R. H. J., R. de Borst, W. a. M. Brekelmans, and M. G. D. Geers.
 * “Gradient-Enhanced Damage Modelling of Concrete Fracture.”
 * Mechanics of Cohesive-Frictional Materials 3, no. 4 (1998): 323–42.
 * https://doi.org/10.1002/(SICI)1099-1484(1998100)3:4<323::AID-CFM51>3.0.CO;2-Z.
 */
class WeaselMaterialPeerlingsConcreteDamage : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  WeaselMaterialPeerlingsConcreteDamage(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// compute scalar isotropic damage based on kappa, softening modulus epsilon_f and an exponential softening law
  std::pair<Real, Vector6r> computeModifiedVonMisesStrain(const Vector6r & strain);

  /// compute scalar isotropic damage based on kappa, softening modulus epsilon_f and an exponential softening law
  std::pair<Real, Real> computeKappa(Real modifiedVonMisesStrain, Real kappaOld);

  /// compute scalar isotropic damage based on kappa, softening modulus epsilon_f and an exponential softening law
  std::pair<Real, Real> computeDamage(Real kappa);

  /// base name for multiple materials on a block
  const std::string _base_name;

  /// Young's modulus
  const Real _E;
  /// Poisson ratio
  const Real _nu;
  /// Poisson ratio
  const Real _kappa_0;
  /// Damage function parameter
  const Real _alpha;
  /// Damage function parameter
  const Real _beta;
  /// Ratio strength conpression tension
  const Real _k;
  /// elastic stiffness
  const Matrix6r _Cel;

  const MaterialProperty<Real> & _kappa_old;
  MaterialProperty<Real> & _kappa;

  MaterialProperty<Real> & _omega;

  const MaterialProperty<RankTwoTensor> & _strain;
  MaterialProperty<RankTwoTensor> & _stress;
  MaterialProperty<RankFourTensor> & _dstress_dstrain;
};
