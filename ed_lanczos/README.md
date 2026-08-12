# ed_lanczos

Exact-diagonalization / Lanczos code for spin-1/2 lattice models (Kitaev–Heisenberg–Γ
type Hamiltonians), with both real-space and momentum-space (k-space, symmetry-reduced)
solvers.

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

If Eigen is installed system-wide (e.g. via `apt install libeigen3-dev` or a
Conan/vcpkg toolchain file), `find_package(Eigen3)` will pick it up
automatically and you can drop `-DEIGEN3_INCLUDE_DIR`.
