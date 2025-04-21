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



//demonstrates a very simple example of numerically integrating (x+1)^4 
//from -2.0 to 3.0 using a Gauss quadrature grid of order 3
#include <cstdio>
#include <cstdlib> 
#include <cstring> 
#include <cmath> 

#include <iostream>
#include <vector>


#include <quadgrid/gauss_legendre_grid.hpp>
using namespace quadgrid;


int main()
{
  const size_t N = 3;

  std::vector<double> x;
  std::vector<double> w;
  double a = -2.0;
  double b = 3.0;


  //1) get grid
  if (!gaussLegendreGrid (N, x, w, a, b))
  {
    std::cout << "gaussLegendreGrid failed for N = " << N << "\n";
    exit(0);
  }

  //2) numerical integral
  double sum1 = 0;
  for (size_t i = 0; i < N; i++)
    sum1 += w[i]*(x[i]+1.0)*(x[i]+1.0)*(x[i]+1.0)*(x[i]+1.0);

  //3) analytical integral
  double sum2 = 1.0/5.0*( (b+1)*(b+1)*(b+1)*(b+1)*(b+1) - 
                          (a+1)*(a+1)*(a+1)*(a+1)*(a+1) );

  char sTemp[500];
  sprintf(sTemp, "integral of (x+1)^4 from a = %.2f to b = %.2f\n", a, b);
  std::cout << sTemp;

  sprintf(sTemp, "sum = %f (analytical)\n", sum2);
  std::cout << sTemp;

  sprintf(sTemp, "sum = %f (numerical)\n", sum1);
  std::cout << sTemp;

  sprintf(sTemp, "error = %.2le (absolute) %.2le (relative)\n", 
    fabs(sum1-sum2), fabs(sum1-sum2)/fabs(sum2));
  std::cout << sTemp;

  return 1;
}





