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
#include <cmath>

#include <complex>
#include <vector>

#include <quadgrid/spherical_harmonic.hpp>
#include <quadgrid/constant.hpp>

namespace quadgrid
{
struct sphereHarmonicTable
{
  size_t lmax = 0;
  double  Y00;
  double  Y10;
  double  Y11;

  std::vector<double> C1;
  std::vector<double> C2;
  std::vector<double> C3;
  std::vector<double> C4;
  std::vector<double> B1;
  std::vector<double> B2;

  void initialize (size_t lmaxInput)
  {
    if ((lmaxInput < lmax) && (lmax != 0))
      return;

    lmax = lmaxInput;

    C1.resize (lmax+1);
    C2.resize (lmax+1);
    C3.resize (lmax+1);
    C4.resize (lmax+1);

    B1.resize (2*lmax+1);
    B2.resize (2*lmax+1);


    for (size_t l = 0; l <= lmax; l++)
    {
      C1[l] = -sqrt( (2.0*l+3)/(2.0*l+2) );
      C2[l] =  sqrt(2.0*l+3);
      C3[l] =  sqrt( (2.0*l+1.0)*(2.0*l+3.0) );
      C4[l] =  sqrt( (2.0*l+3.0)/(2.0*l-1.0) );
    }

    for (size_t l = 0; l <= 2*lmax; l++)
    {
      B1[l] = sqrt( (double) l);
      B2[l] = sqrt( (double) l/(l+1.0) );
    }


    Y00 = sqrt(1.0/4.0/Pi);
    Y10 = Y00*C2[0];
    Y11 = Y00*C1[0];
  }

};




sphereHarmonicTable& getSphereHarmonicTable ()
{
  static sphereHarmonicTable table;
  return table;
}



void sphereHarmonic (std::vector<std::complex<double>>& Ylm, const size_t lmax, 
  const double rUnit[3])
{
  auto& table = getSphereHarmonicTable();
  table.initialize(lmax);


  const size_t sizeYlm = sphereHarmonicArraySize (lmax);
  Ylm.resize(sizeYlm);


  Ylm[0] = table.Y00;
  if (lmax == 0)
    return;


  std::complex<double> r01 = std::complex<double>(rUnit[0], rUnit[1]);

  Ylm[1] = table.Y10*rUnit[2];
  Ylm[2] = table.Y11*r01;

  size_t index1, index2, index3;

  size_t l_index    = 3;
  size_t l_index_m1 = 1;
  size_t l_index_m2 = 0;

  for (size_t l = 2; l <= lmax; l++)
  {
    for (size_t m = 0; m <= l-2; m++)
    {
      index1 = l_index    + m;
      index2 = l_index_m1 + m;
      index3 = l_index_m2 + m;

      double fact1  = table.C3[l-1]/table.B1[l+m]/table.B1[l-m]*rUnit[2];
      double fact2  = table.C4[l-1]*table.B2[l-1+m]*table.B2[l-1-m];
      Ylm[index1]   = fact1*Ylm[index2] - fact2*Ylm[index3];
    }

    //Downward Recursion
    index1      = l_index    + l-1;
    index2      = l_index_m1 + l-1;
    Ylm[index1] = table.C2[l-1]*rUnit[2]*Ylm[index2];

    index1      = l_index    + l;
    Ylm[index1] = table.C1[l-1]*r01*Ylm[index2];

    l_index     += l + 1;
    l_index_m1  += l;
    l_index_m2  += l - 1;
  }
}

void sphereHarmonic (std::vector<std::complex<double>>& Ylm, const size_t lmax, 
  const double theta, const double phi)
{
  double rUnit[3];
  
  const double cos_theta = cos(theta);
  const double sin_theta = sin(theta);
  const double cos_phi   = cos(phi);
  const double sin_phi   = sin(phi);

  rUnit[0] = sin_theta*cos_phi;
  rUnit[1] = sin_theta*sin_phi;
  rUnit[2] = cos_theta;

  return sphereHarmonic (Ylm, lmax, rUnit);
}







}//end namespace quadgrid



