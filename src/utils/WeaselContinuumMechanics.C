#include "WeaselContinuumMechanics.h"

namespace Weasel
{

Vector6r
voigtFromRankTwoTensor(const RankTwoTensor & tensor, bool multiplyShearTermsX2)
{
  Vector6r voigt;

  voigt(0) = tensor(0, 0);
  voigt(1) = tensor(1, 1);
  voigt(2) = tensor(2, 2);
  voigt(3) = tensor(0, 1);
  voigt(4) = tensor(0, 2);
  voigt(5) = tensor(1, 2);

  if (multiplyShearTermsX2)
  {
    voigt(3) *= 2;
    voigt(4) *= 2;
    voigt(5) *= 2;
  }

  return voigt;
}

RankTwoTensor
rankTwoTensorFromVoigt(const Vector6r & r2tInVoigt, bool divideShearTermsBy2)
{
  const auto & v = r2tInVoigt;
  RankTwoTensor r2t = RankTwoTensor(v(0),
                                    v(1),
                                    v(2),
                                    divideShearTermsBy2 ? v(5) / 2 : v(5),
                                    divideShearTermsBy2 ? v(4) / 2 : v(4),
                                    divideShearTermsBy2 ? v(3) / 2 : v(3));

  return r2t;
}

RankFourTensor
rankFourTensorFromVoigt(const Matrix6r & r4tInVoigt)
{
  const static std::array<std::array<unsigned int, 3>, 3> comp2vgt{
      {{0, 3, 4}, {3, 1, 5}, {4, 5, 2}}};

  RankFourTensor r4t;

  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
      for (unsigned k = 0; k < 3; ++k)
        for (unsigned l = 0; l < 3; ++l)
        {
          r4t(i, j, k, l) = r4tInVoigt(comp2vgt[i][j], comp2vgt[k][l]);

          /* if ( i != j && _divide_shear_terms_by_2_ij ) */
          /*   _the_rank_four_tensor[_qp]( i, j, k, l ) *= 0.5; */
          /* if ( k != l && _divide_shear_terms_by_2_kl ) */
          /*   _the_rank_four_tensor[_qp]( i, j, k, l ) *= 2; */
        }

  return r4t;
}

namespace Elasticity
{

Matrix6r
stiffnessTensor(const Real E, const Real nu)
{
  Matrix6r C;
  // clang-format off
            C <<  (1-nu),     nu,     nu,          0,          0,          0,
                      nu, (1-nu),     nu,          0,          0,          0,
                      nu,     nu, (1-nu),          0,          0,          0,
                      0,       0,      0, (1-2*nu)/2,          0,          0,
                      0,       0,      0,          0, (1-2*nu)/2,          0,
                      0,       0,      0,          0,          0, (1-2*nu)/2;
  // clang-format on
  C *= E / ((1 + nu) * (1 - 2 * nu));
  return C;
}

Matrix6r
complianceTensor(const Real E, const Real nu)
{
  Matrix6r CInv;
  const Real G = E / (2 * (1 + nu));
  // clang-format off
              CInv <<   1./E,  -nu/E, -nu/E,     0,    0,    0,
                       -nu/E,   1./E, -nu/E,     0,    0,    0,
                       -nu/E,  -nu/E,  1./E,     0,    0,    0,
                           0,      0,     0,  1./G,    0,    0,
                           0,      0,     0,     0, 1./G,    0,
                           0,      0,     0,     0,    0, 1./G;
  // clang-format on
  return CInv;
}

}

namespace Invariants
{
const Vector6r P = (Vector6r() << 1, 1, 1, 2, 2, 2).finished();
const Vector6r PInv = (Vector6r() << 1, 1, 1, .5, .5, .5).finished();

const Vector6r I = (Vector6r() << 1, 1, 1, 0, 0, 0).finished();
const Vector6r IHyd = (Vector6r() << 1. / 3, 1. / 3, 1. / 3, 0, 0, 0).finished();

const Matrix6r IDev = (Matrix6r() <<
                           // clang-format off
        2./3,    -1./3,   -1./3,    0,  0,  0,
        -1./3,   2./3,    -1./3,    0,  0,  0,
        -1./3,   -1./3,   2./3,     0,  0,  0,
        0,          0,      0,      1,  0,  0,
        0,          0,      0,      0,  1,  0,
        0,          0,      0,      0,  0,  1).finished();
// clang-format on

Real
I1(const Vector6r & stress)
{
  return stress.head(3).sum();
}

Real
I2(const Vector6r & stress)
{
  const Vector6r & s = stress;

  return s(0) * s(1) + s(1) * s(2) + s(2) * s(0) - s(3) * s(3) - s(4) * s(4) - s(5) * s(5);
}
Real
I3(const Vector6r & stress)
{
  const Vector6r & s = stress;
  return s(0) * s(1) * s(2) + 2 * s(3) * s(4) * s(5) - s(0) * s(5) * s(5) - s(1) * s(4) * s(4) -
         s(2) * s(3) * s(3);
}

Real
J2(const Vector6r & stress)
{

  Real I1_ = I1(stress);
  Real I2_ = I2(stress);
  Real res = (1. / 3) * I1_ * I1_ - I2_;
  return std::max(res, 0.0);
}
Real
J3(const Vector6r & stress)
{
  Real I1_ = I1(stress);
  Real I2_ = I2(stress);
  Real I3_ = I3(stress);

  return (2. / 27) * pow(I1_, 3) - (1. / 3) * I1_ * I2_ + I3_;
}

Vector6r
dJ2_dStress(const Vector6r & stress)
{
  return P.array() * (IDev * stress).array();
}

}

}
