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

#ifndef QUADGRID_UNIT_SPHERE_GRID_GAUSS_LEGENDRE_HPP
#define QUADGRID_UNIT_SPHERE_GRID_GAUSS_LEGENDRE_HPP


#include <vector>
#include <cstddef>

namespace quadgrid
{

/// \brief Generates a spherical quadrature grid using Gauss-Legendre quadrature in the polar direction.
/// \param lmax Approximate order of the grid (max degree of spherical harmonics to integrate accurately).
/// \param N Output number of polar grid points (will be approximately `lmax / 2`).
/// \param theta Output vector of polar angles \( \theta \in [0, \pi] \) (length N).
/// \param w_theta Output vector of associated polar weights (length N). Azimuthal weights are folded into w_theta.
/// \param phi Output vector of azimuthal angles \( \phi \in [0, 2\pi] \) (length 2×N).
/// \return `true` on success; `false` if `lmax` is out of bounds (e.g., too large or unsupported).
/// \note This grid uses Gauss-Legendre sampling in the polar direction and uniform sampling in the azimuthal direction.
/// \note `lmax` should not exceed 2000 for accurate results.
bool unitSphereGaussLegendre(const size_t lmax,
                             size_t& N,
                             std::vector<double>& theta,
                             std::vector<double>& w_theta,
                             std::vector<double>& phi);

} // namespace quadgrid




#endif //QUADGRID_UNIT_SPHERE_GRID_GAUSS_LEGENDRE_HPP



