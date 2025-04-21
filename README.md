# quadgrid

**quadgrid** is a lightweight C++ library providing high-accuracy quadrature grids for scientific integration over the interval and the unit sphere. 

---

## Features

- Precomputed **Gauss-Legendre** quadrature grids (1D) for N = 1-1000
- Full set of **Lebedev** unit sphere grids (for spherical integration)
- Custom **spherical Gauss-Legendre** grid (latitudinal and longitudinal sampling)
- Supporting utilities: Legendre polynomials and real/complex spherical harmonics for testing and convergence analysis
- Header-only interface with minimal dependencies
- Numerically verified: spherical harmonics integration errors ≤ **3e-14**

---

## Build Instructions

mkdir build && cd build
cmake ..
make


## Example
#include <quadgrid/gauss_legendre_grid.hpp>

std::vector<double> x, w;
quadgrid::gaussLegendreGrid (64, x, w, 2.0, 5.0);
// creates a Gauss Legendre grid of order 64 for the [2.0, 5.0] interval


📖 [View the Documentation](https://elkingd.github.io/quadgrid/index.html)


