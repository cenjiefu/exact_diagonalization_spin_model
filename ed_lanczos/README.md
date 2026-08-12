# About

Exact-diagonalization / Lanczos code for spin-1/2 lattice models, with both real-space and momentum-space (k-space, symmetry-reduced)
solvers.
Can define arbitrary lattice with any spin interactions between sites. The example is the extended Kitaev models on honeycomb lattices of 6,12,24 sites.

## Dependencies

All are required to build:

- A C++17 compiler
- [Eigen](https://eigen.tuxfamily.org/) 3.4+ (dense/sparse linear algebra)
- [EigenRand](https://github.com/bab2min/EigenRand) (header-only)
- [ezarpack](https://github.com/krivenko/ezarpack) (header-only wrapper around ARPACK)
- ARPACK (e.g. [arpack-ng](https://github.com/opencollab/arpack-ng)) + BLAS/LAPACK
- OpenMP

Eigen, EigenRand, and ezarpack are header-only and are **not vendored in this
repo**. If they aren't already visible to your compiler's default include
path, point CMake at them

If Eigen is installed system-wide (e.g. with "module load arpack-ng Eigen3" in HPC clusters), `find_package(Eigen3)` will pick it up
automatically

## Usage
change line 58 and 59 in CMakeLists.txt to your local path to eigen and ezarpack.
then build with cmake
