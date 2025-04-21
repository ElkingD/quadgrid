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




//tests the orthogonality of the Legendre polynomials by numerical integration
//with the Gauss-Legendre grids


#include <cstdio>
#include <cstdlib> 
#include <cstring> 
#include <cmath> 

#include <iostream>
#include <vector>

#include <quadgrid/legendre.hpp>
#include <quadgrid/gauss_legendre_grid.hpp>
using namespace quadgrid;


int main()
{
  //set up grid size arrays (1, 2, .. 100) and (110, 120, .. 1000)
  std::vector<size_t> arrayOrder;
  for (size_t N = 1; N <= 100; N++)
    arrayOrder.push_back(N);
  for (size_t N = 110; N <= 1000; N += 10)
    arrayOrder.push_back(N);

  
  for (size_t i = 0; i < arrayOrder.size(); i++)
  {
    const size_t N = arrayOrder[i];

    //1) get grid coordinates (x) and weights (w)
    std::vector<double> x;
    std::vector<double> w;
    if (!gaussLegendreGrid (N, x, w, -1.0, 1.0))
    {
      std::cout << "gaussLegendreGrid failed for N = " << N << "\n";
      exit(0);
    }

     
    //2) initialize numerical integral sum to zero for each order n = 0,1,..N
    std::vector<double> P(N+1);
    std::vector<double> sum(N+1);
    for (size_t i = 0; i <= N; i++)
      sum[i] = 0;


    //3) calculate numerical integral for each order n = 0,1,..N
    for (size_t i = 0; i < N; i++)
    {
      legendrePoly (P, x[i], N);
      for (size_t n = 0; n <= N; n++)
        sum[n] += w[i]*P[n];
    }

    //4) normalize by (2n+1)/2 and measure errors.
    //sum[n] = 1 (for n = 0) and 0 (for n != 0)

    double maxError = 0.0;
    double rmsError = 0.0;
    for (size_t n = 0; n <= N; n++)
    {
      sum[n] *= (2.0*n+1.0)/2.0;

      double error = fabs(sum[n]);
      if (n == 0)
        error = fabs(sum[n]-1.0);

      if (maxError < error)
        maxError = error;
      rmsError += error*error;
    }
    rmsError /= (N+1);
    rmsError = sqrt(rmsError);



    {
      char sTmp[500];
      sprintf(sTmp, "Grid %4lu maxError = %.2le  rmsError = %.2le ((2*m+1)/2*<1|Pm> = delta(m,0))\n", 
        N, maxError, rmsError);
      std::cout << sTmp;
    }
  }




  return 1;
}





