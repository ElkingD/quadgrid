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

#include <vector>

#include <quadgrid/legendre.hpp>


namespace quadgrid
{
void legendrePoly (std::vector<double>& P, const double x, const size_t N)
//calculates P[N+1] up to N
{
  P.resize(N+1);

  if (N == 0) 
  {
    P[0] = 1.0;
    return;
  }
  else if (N == 1)
  {
    P[0] = 1.0;
    P[1] = x;
    return;
  }

  P[0] = 1.0;
  P[1] = x;

  const double x_2 = 2.0*x;
  for (size_t n = 1; n < N; n++)
    P[n+1] = (x_2*P[n] - P[n-1]) - (x*P[n] - P[n-1])/(n+1.0);

}
void legendrePoly (std::vector<double>& P, std::vector<double>& dPdx, 
  const double x, const size_t N)
//calculates P[N+1] and derivative dPdx[N+1] up to N
{
  P.resize(N+1);
  dPdx.resize(N+1);

  if (N == 0)
  {
    P[0]    = 1.0;
    dPdx[0] = 0.0;
    return;
  }
  else if (N == 1)
  {
    P[0]    = 1.0;
    P[1]    = x;
    dPdx[0] = 0.0;
    dPdx[1] = 1.0;
    return;
  }

  P[0]    = 1.0;
  P[1]    = x;
  const double x_2 = 2.0*x;
  for (size_t n = 1; n < N; n++)
    P[n+1] = (x_2*P[n] - P[n-1]) - (x*P[n] - P[n-1])/(n+1.0);

  dPdx[0] = 0.0;
  dPdx[1] = 1.0;
  for (size_t n = 2; n <= N; n++)
    dPdx[n] = n*P[n-1] + x*dPdx[n-1];
}


}//end namespace quadgrid


