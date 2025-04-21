// SPDX-License-Identifier: Apache-2.0
//
// quadgrid - High-accuracy quadrature grids for scientific computing
// Copyright 2025 Denny Elking
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//

#ifndef QUADGRID_LEGENDRE_HPP
#define QUADGRID_LEGENDRE_HPP

/// \file
/// \brief Functions for computing Legendre polynomials.
#include <vector>
#include <cstddef>

namespace quadgrid
{

/// \brief Computes Legendre polynomials \( P_0(x), P_1(x), \dots, P_N(x) \).
/// \param P Output vector that will contain \( N+1 \) values: P[0] through P[N].
/// \param x Input value at which to evaluate the polynomials.
/// \param N Maximum order of the Legendre polynomial.
/// \note The output vector P is resized if necessary.
void legendrePoly(std::vector<double>& P, const double x, const size_t N);

/// \brief Computes Legendre polynomials and their first derivatives at point \( x \).
/// \param P Output vector for Legendre polynomials \( P_0(x), \dots, P_N(x) \) (size N+1).
/// \param dPdx Output vector for first derivatives \( \frac{d}{dx}P_0(x), \dots, \frac{d}{dx}P_N(x) \).
/// \param x Input value at which to evaluate the polynomials and derivatives.
/// \param N Maximum order of the Legendre polynomial.
/// \note Both output vectors P and dPdx are resized if necessary.
void legendrePoly(std::vector<double>& P,
                  std::vector<double>& dPdx,
                  const double x,
                  const size_t N);



}//end namespace quadgrid




#endif //QUADGRID_LEGENDRE_HPP


