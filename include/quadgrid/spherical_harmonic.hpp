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

#ifndef QUADGRID_SPHERICAL_HPPARMONIC_HPP
#define QUADGRID_SPHERICAL_HPPARMONIC_HPP




#ifndef QUADGRID_SPHERICAL_HARMONIC_HPP
#define QUADGRID_SPHERICAL_HARMONIC_HPP

#include <vector>
#include <complex>
#include <cstddef>

namespace quadgrid
{

/// \brief Returns the number of complex spherical harmonic values Y(l,m), using a compact 1D layout.
/// \param lmax The maximum degree of the spherical harmonic.
/// \return Total number of Y(l,m) values for l >= 0 and l >= m >= 0.
/// \note Array stores only non-negative m; negative m values are recovered using symmetry.
inline size_t sphereHarmonicArraySize(size_t lmax)
{
  return ((lmax + 1) * (lmax + 2)) / 2;
}

/// \brief Computes the flat array index for Y(l,m) where l >= m >= 0.
/// \param l Degree index.
/// \param m Order index.
/// \return Linear index for accessing Y(l,m) in a flattened array.
/// \note Assumes a packed layout: [Y(0,0),Y(1,0),Y(1,1),Y(2,0),Y(2,1) ...]
inline size_t sphereHarmonicArrayIndex(size_t l, size_t m)
{
  return (l * (l + 1)) / 2 + m;
}

/// \brief Computes complex spherical harmonics Y(l,m) up to degree lmax using a Cartesian unit vector.
/// \param Ylm Output vector containing Y(l,m) values in compact layout (size = sphereHarmonicArraySize(lmax)).
/// \param lmax Maximum degree to compute.
/// \param rUnit Unit vector [x, y, z] on the sphere.
/// \note Computes only Y(l,m) for l >= m >= 0
void sphereHarmonic(std::vector<std::complex<double>>& Ylm,
                    const size_t lmax,
                    const double rUnit[3]);

/// \brief Computes complex spherical harmonics Y(l,m) up to degree lmax using spherical coordinates.
/// \param Ylm Output vector containing Y(l,m) values in compact layout (size = sphereHarmonicArraySize(lmax)).
/// \param lmax Maximum degree to compute.
/// \param theta Polar angle 0 <= theta <= Pi (from +z axis).
/// \param phi Azimuthal angle 0 <= phi <= 2*Pi (around z axis).
/// \note Computes only Y(l,m) for l >= m >= 0
void sphereHarmonic(std::vector<std::complex<double>>& Ylm,
                    const size_t lmax,
                    const double theta,
                    const double phi);

} // namespace quadgrid

#endif // QUADGRID_SPHERICAL_HARMONIC_HPP





#endif //QUADGRID_SPHERICAL_HPPARMONIC_HPP












