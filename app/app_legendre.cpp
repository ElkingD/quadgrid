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

#include <cstdio>
#include <cstdlib> 
#include <cstring> 
#include <cmath> 

#include <quadgrid/legendre.hpp>


using namespace quadgrid;


int main()
{
  double x = 0.3;
  const size_t N = 10;


  std::vector<double> P;
  std::vector<double> dPdx;


  legendrePoly (P, x, N);
  //calculates legendre polynomials P[N+1] up to N

  for (size_t n = 0; n <= N; n++)
    printf("P[%3lu](%6.3f) = %10.6f\n", n, x, P[n]);

  legendrePoly (P, dPdx, x, N);

  //calculates legendre polynomials and its first derivative 
  //P[N+1] and dPdx[N+1] up to N


  for (size_t n = 0; n <= N; n++)
    printf("P[%3lu](%6.3f) = %10.6f dPdx[%3lu](%6.3f) = %10.6f\n", 
      n, x, P[n], n, x, dPdx[n]);


  return 1;
}





