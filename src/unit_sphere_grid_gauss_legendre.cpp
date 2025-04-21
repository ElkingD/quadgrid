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

#include <cstring>
#include <cstdlib>
#include <cmath>

#include <iostream>
#include <vector>

#include <quadgrid/unit_sphere_grid_gauss_legendre.hpp>
#include <quadgrid/gauss_legendre_grid.hpp>
#include <quadgrid/constant.hpp>


namespace quadgrid
{


bool unitSphereGaussLegendre (const size_t lmax, size_t& N,
  std::vector<double>& theta, std::vector<double>& w_theta, 
  std::vector<double>& phi)
//input:  lmax
//output: N ~ lmax/2
//        theta[N]
//        w_theta[N]   phi weight is absorbed into w_theta
//        phi[2*N]
{
  if (lmax%2 == 0)
    N = (lmax+2)/2;
  else
    N = (lmax+1)/2;

  //the Gauss Legendre grids are defined up to order 1000
  if (N > 1000)
  {
    std::cout << "Error in unitSphereGaussLegendre\n";
    std::cout << "  lmax = " << lmax << " N = " << N << " > 1000\n";
    return false;
  }

  if (N > 100)
  //for N > 100, Legendre grids only every += 10 up to N == 1000
  {
    if (N%10 != 0)
      N = 10*(N/10 + 1);
  }

  theta.resize(N);
  w_theta.resize(N);
  phi.resize(2*N);


  //weights for u = cos(theta), du = -sin(theta)dTheta 
  // 0 <= theta <= Pi
  //-1 <= u     <= 1
  if (!gaussLegendreGrid (N, theta, w_theta, -1.0, 1.0))
  {
    std::cout << "Error in unitSphereGaussLegendre\n";
    std::cout << "  gaussLegendreGrid failed for N = " << N << "\n";
    return false;
  }

  for (size_t i = 0; i < N; i++)
  {
    theta[i]    = acos(theta[i]);
    w_theta[i] *= Pi/N;
  }

  for (size_t i = 0; i < 2*N; i++)
    phi[i] = Pi/N*i;


  return true;
}

}//end namespace quadgrid



