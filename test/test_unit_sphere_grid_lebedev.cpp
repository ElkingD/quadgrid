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




//for each unit sphere lebedev grid with order lmax, the following test
//numericall integrates the spherical harmonics Y(l,m) for l <= lmax
//due to the orthonormality of Y(l,m), only the integral of Y(0,0) is non-zero


#include <cstdio>
#include <cstdlib> 
#include <cstring> 
#include <cmath> 

#include <iostream>
#include <vector>
#include <complex>

#include <quadgrid/unit_sphere_grid_lebedev.hpp>
#include <quadgrid/spherical_harmonic.hpp>
#include <quadgrid/constant.hpp>
using namespace quadgrid;


int main()
{
  bool PRINT_ALL = 0;

  for (size_t i = 0; i < unitSphereLebedevNumGrid; i++)
  //there are 11 grids.  test each one.
  {
    //1) get grid
    size_t lmax, nPoint;
    std::vector<double> coord;
    std::vector<double> weight;
    if (!unitSphereLebedev (i, lmax, nPoint, coord, weight))
    {
      std::cout << "Error.  unitSphereLebedev failed for " << i << "\n";
      exit(0);
    }

     //2) initial numerical sum to zero for each (l,m)
    const size_t sizeYlm = sphereHarmonicArraySize (lmax);
    std::vector< std::complex<double> > sum(sizeYlm);
    for (size_t n = 0; n < sizeYlm; n++)
      sum[n] = 0;


    //3) calculate spherical harmonics at each point and add to numerical sum
    std::vector< std::complex<double> > Ylm;
    for (size_t n1 = 0; n1 < nPoint; n1++)
    {
      const double *rUnit = &coord[3*n1];
      sphereHarmonic (Ylm, lmax, rUnit);

      for (size_t n2 = 0; n2 < sizeYlm; n2++)
        sum[n2] += weight[n1]*Ylm[n2];
    }    

    //4) normalize by Y00
    const double Y00 = sqrt(1.0/4.0/Pi);
    for (size_t n2 = 0; n2 < sizeYlm; n2++)
      sum[n2] *= Y00;
    

    //5) only sum[0] == 1.0 and the other's should be zero
    double maxError = 0;
    double rmsError = 0;
    for (size_t n = 0; n < sizeYlm; n++)
    {
      if (PRINT_ALL)
      {
        char sTmp[500];
        sprintf(sTmp, "sum[%4lu] = %15.12f + %15.12fi\n", 
          n, sum[n].real(), sum[n].imag());
        std::cout << sTmp;
      }

      double error = abs(sum[n]);
      if (n == 0)
        error = abs(sum[n]-1.0);


      if (maxError < error)
        maxError = error;
      rmsError += error*error;
    }
    rmsError /= sizeYlm;
    rmsError = sqrt(rmsError);

    {
      char sTmp[500];
      sprintf(sTmp, "grid = %4lu lmax = %4lu maxError = %.2le rmsError = %.2le\n",
        i, lmax, maxError, rmsError);
      std::cout << sTmp;
    }
  }







  return 1;
}





