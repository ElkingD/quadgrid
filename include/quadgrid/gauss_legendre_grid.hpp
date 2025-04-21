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


#ifndef QUADGRID_GAUSS_LEGENDRE_GRID_HPP
#define QUADGRID_GAUSS_LEGENDRE_GRID_HPP


#include <vector>
#include <cstddef>

namespace quadgrid
{
/// \brief Computes Gauss-Legendre quadrature nodes and weights.
/// \param N The number of quadrature points (order).
///          Must be one of: {1, 2, 3.. 100} or {110, 120, 130, .. 1000}
/// \param x Output vector to store the quadrature nodes (size N).
/// \param w Output vector to store the corresponding weights (size N).
/// \param a Lower bound of integration interval [a, b]
/// \param b Upper bound of integration interval [a, b]
/// \note The weights and nodes are for integration over [a, b].
bool gaussLegendreGrid (const size_t N, std::vector<double>& x, 
  std::vector<double>& w, const double a, const double b);
//input:  N = order, [a, b] interval
//output: x[N] and w[N] = coordinates and weights
//        Integral{ f(x) } over [a,b] = Sum{ f(x[i])*w[i] } from i = 0 to N - 1


}//end namespace quadgrid




#endif //QUADGRID_GAUSS_LEGENDRE_GRID_HPP


