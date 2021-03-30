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
#include "WeaselMaterialPeerlingsConcreteDamage.h"

class WeaselMaterialPeerlingsConcreteGradientDamage : public WeaselMaterialPeerlingsConcreteDamage
{
public:
  static InputParameters validParams();

  WeaselMaterialPeerlingsConcreteGradientDamage(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  const Real & _the_nonlocal_radius;
  const VariableValue & _equivalent_strain_nonlocal;

  MaterialProperty<Real> & _equivalent_strain_local;
  MaterialProperty<Real> & _dequivalent_strain_local_dequivalent_strain_nonlocal;
  MaterialProperty<RankTwoTensor> & _dequivalent_strain_local_dstrain;
  MaterialProperty<RankTwoTensor> & _dstress_dequivalent_strain_nonlocal;

  MaterialProperty<Real> & _nonlocal_radius;
};
