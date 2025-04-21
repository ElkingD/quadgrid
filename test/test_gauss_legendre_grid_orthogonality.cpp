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
  const size_t MaxOrder = 1000;

  bool PRINT_ALL = 0;

  std::vector<double> x;
  std::vector<double> w;
  std::vector<double> Ptmp;


  std::vector<double> Parray((MaxOrder+1)*MaxOrder);

  //set up grid size arrays (1, 2, .. 100) and (110, 120, .. 1000)
  std::vector<size_t> arrayOrder;
  for (size_t N = 1; N <= 100; N++)
    arrayOrder.push_back(N);
  for (size_t N = 110; N <= 1000; N += 10)
    arrayOrder.push_back(N);

  
  for (size_t i = 0; i < arrayOrder.size(); i++)
  {
    const size_t N = arrayOrder[i];

    //1) get grid
    if (!gaussLegendreGrid (N, x, w, -1.0, 1.0))
    {
      std::cout << "gaussLegendreGrid failed for N = " << N << "\n";
      exit(0);
    }

    if (PRINT_ALL)
      std::cout << "Grid " << N << "\n";

      
    //2) calculate legendre polynomials for each grid point and store
    for (size_t i = 0; i < N; i++)
    {
      legendrePoly (Ptmp, x[i], N);
      for (size_t n = 0; n <= N; n++)
      {
        const size_t index = n*N+i;
        Parray[index]   = Ptmp[n];
      }
    }

    //3) test orthonormality of Legendre polynomials (multiplied by (2*n+1)/2)
    double maxError = 0.0;
    double rmsError = 0.0;
    int count = 0;
    for (size_t n = 1; n <= N; n++)
    {
      for (size_t m = 1; m <= N-n; m++)
      {
        double fnume = 0.0;
        const double *Pn = &Parray[n*N];
        const double *Pm = &Parray[m*N];
        for (size_t i = 0; i < N; i++)
        {
          fnume += w[i]*Pn[i]*Pm[i];
        }
        fnume *= (2.0*n+1.0)/2;
        double fanal = 0.0;
        if (n == m)
          fanal = 1.0;

        double error = fabs(fnume-fanal);
        if (maxError < error) maxError = error;
        rmsError += error*error;
        count++;

        if (error > 1.0E-10)
          PRINT_ALL = 1;

        if (PRINT_ALL)
        {
          char sTmp[500];
          const double C = 2.0/(2.0*n+1.0);

          sprintf(sTmp, "  <P%-3lu|P%-3lu> = %20.15f (anal)\n", n, m, C*fanal);
          std::cout << sTmp;

          sprintf(sTmp, "  <P%-3lu|P%-3lu> = %20.15f (nume)\n", n, m, C*fnume);
          std::cout << sTmp;

          sprintf(sTmp, "  error = %.2le\n\n", C*error);
          std::cout << sTmp;
        }

        if (error > 1.0E-10)
        {
          char sTmp[500];
          sprintf(sTmp, "Error. error = %.2le > 1.0E-10\n", error);
          std::cout << sTmp;
          exit(0);
        }

      }
    }
    rmsError /= count;
    rmsError  = sqrt(rmsError);

    {
      char sTmp[500];
      sprintf(sTmp, "Grid %4lu maxError = %.2le  rmsError = %.2le (all <Pn|Pm> nume integrals for n+m <= %4lu)\n", 
        N, maxError, rmsError, N);
      std::cout << sTmp;
    }
  }




  return 1;
}





