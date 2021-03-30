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
 * Computes the divergence of 2.order stress tensors
 */

class WeaselStressDivergence : public DerivativeMaterialInterface<Kernel>
{
public:
  static InputParameters validParams();

  WeaselStressDivergence(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  Real computeQpJacobianDisplacement(unsigned int comp_i, unsigned int comp_j);
  Real computeQpJacobianNonlocalDamage(unsigned int comp_i);

  /// Base name of the material system that this kernel applies to
  const std::string _base_name;

  /// The tensor
  const MaterialProperty<RankTwoTensor> & _stress;

  /// Derivatives of the w.r.t. strain increment
  const MaterialProperty<RankFourTensor> & _dstress_dstrain;

  /// Is there a nonlocal damage field?
  const bool _has_nonlocal_damage;

  /// Derivatives of the w.r.t. the nonlocal damage driving field
  const MaterialProperty<RankTwoTensor> & _dstress_dk;

  /// An integer corresponding to the direction this kernel acts in
  const unsigned int _component;

  /// Coupled displacement variables
  unsigned int _ndisp;
  /// Displacement variables IDs
  std::vector<unsigned int> _disp_var;

  /// The MOOSE variable number of the nonlocal damage variable
  unsigned int _nonlocal_damage_var;
};
