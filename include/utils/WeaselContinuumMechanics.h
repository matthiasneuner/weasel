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
#include "Eigen/Core"
#include "Moose.h"
#include "MooseTypes.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

using Vector6r = Eigen::Matrix<Real, 6, 1>;
using Matrix6r = Eigen::Matrix<Real, 6, 6>;

using VectorXr = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
using MatrixXr = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

namespace Weasel
{

Vector6r voigtFromRankTwoTensor(const RankTwoTensor & tensor, bool multiplyShearTermsX2);

RankTwoTensor rankTwoTensorFromVoigt(const Vector6r & r2tInVoigt, bool divideShearTermsBy2);

RankFourTensor rankFourTensorFromVoigt(const Matrix6r & r4tInVoigt);

namespace Elasticity
{

Matrix6r stiffnessTensor(const Real E, const Real nu);

Matrix6r complianceTensor(const Real E, const Real nu);

}

namespace Invariants
{

Real I1(const Vector6r & stress);
Real I2(const Vector6r & stress);
Real I3(const Vector6r & stress);

Real J2(const Vector6r & stress);
Real J3(const Vector6r & stress);

Vector6r dJ2_dStress(const Vector6r & stress);

}

}
