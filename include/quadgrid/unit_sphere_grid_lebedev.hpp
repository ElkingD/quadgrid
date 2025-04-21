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

#ifndef QUADGRID_UNIT_SPHERE_GRID_LEBEDEV_HPP
#define QUADGRID_UNIT_SPHERE_GRID_LEBEDEV_HPP


#include <vector>

namespace quadgrid
{
// \brief Number of supported Lebedev grids.
const size_t unitSphereLebedevNumGrid = 11;

/// \brief Maximum spherical harmonic order (lmax) supported for each grid.
const size_t unitSphereLebedevLmax[unitSphereLebedevNumGrid] =
{
  9, 11, 17, 23, 25, 27, 29, 41, 47, 53, 131
};

/// \brief Number of integration points on the sphere for each Lebedev grid.
const size_t unitSphereLebedevNumPoint[unitSphereLebedevNumGrid] =
{
  38, 50, 110, 194, 230, 266, 302, 590, 770, 974, 5810
};

/// \brief Retrieves a Lebedev unit sphere grid by index.
/// \param index Index from 0 to 10 identifying the grid size and order.
/// \param lmax Output: maximum spherical harmonic degree integrated accurately by this grid.
/// \param nPoint Output: number of quadrature points on the unit sphere.
/// \param coord Output: flattened coordinate array of size nPoint × 3 (x, y, z for each point).
/// \param weight Output: quadrature weights for each point (size nPoint).
/// \return `true` if the index is valid; `false` if out of range.
/// \note The integration error has been benchmarked by integrating spherical harmonics
///       \( Y_\ell^m \) over the sphere. Only the \( Y_0^0 \) integral should be nonzero.
///       Errors for each grid index are listed below.
///
/// \verbatim
///   grid =  0  lmax =   9   maxError = 2.58e-16  rmsError = 6.57e-17
///   grid =  1  lmax =  11   maxError = 6.07e-16  rmsError = 1.20e-16
///   grid =  2  lmax =  17   maxError = 2.22e-15  rmsError = 2.97e-16
///   grid =  3  lmax =  23   maxError = 2.66e-15  rmsError = 2.52e-16
///   grid =  4  lmax =  25   maxError = 9.14e-15  rmsError = 8.19e-16
///   grid =  5  lmax =  27   maxError = 2.87e-15  rmsError = 4.73e-16
///   grid =  6  lmax =  29   maxError = 2.78e-15  rmsError = 2.36e-16
///   grid =  7  lmax =  41   maxError = 3.11e-15  rmsError = 2.14e-16
///   grid =  8  lmax =  47   maxError = 6.44e-15  rmsError = 2.33e-16
///   grid =  9  lmax =  53   maxError = 1.10e-14  rmsError = 3.39e-16
///   grid = 10  lmax = 131   maxError = 2.73e-14  rmsError = 3.21e-16
/// \endverbatim
bool unitSphereLebedev(const size_t index,
                       size_t& lmax,
                       size_t& nPoint,
                       std::vector<double>& coord,
                       std::vector<double>& weight);




}//end namespace quadgrid




#endif //QUADGRID_UNIT_SPHERE_GRID_LEBEDEV_HPP


