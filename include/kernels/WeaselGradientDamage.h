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

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"

/**
 * Computes the Helmholtz like equation for gradient-enhanced damage formulations
 */

class WeaselGradientDamage : public DerivativeMaterialInterface<Kernel>
{
public:
  static InputParameters validParams();

  WeaselGradientDamage(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  Real computeQpJacobianDisplacement(unsigned int comp_u);
  Real computeQpJacobianNonlocalDamage();

  /// Base name of the material system that this kernel applies to
  const std::string _base_name;

  const MaterialProperty<Real> & _local_damage;
  const MaterialProperty<Real> & _nonlocal_radius;

  /// Derivatives of the w.r.t. strain increment
  const MaterialProperty<RankTwoTensor> & _dlocal_damage_dstrain;

  /// Derivatives of the w.r.t. nonlocal damage
  const MaterialProperty<Real> & _dlocal_damage_dnonlocal_damage;

  /// Coupled displacement variables
  unsigned int _ndisp;
  /// Displacement variables IDs
  std::vector<unsigned int> _disp_var;

  /// The MOOSE variable number of the nonlocal damage variable
  unsigned int _nonlocal_damage_var;
};
