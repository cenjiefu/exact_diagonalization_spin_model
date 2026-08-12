# ed_lanczos

Exact-diagonalization / Lanczos code for spin-1/2 lattice models (Kitaev–Heisenberg–Γ
type Hamiltonians), with both real-space and momentum-space (k-space, symmetry-reduced)
solvers.

This is a reorganization of a working research codebase into a standard C++
library/app layout. **No algorithmic logic was changed** — only file layout,
include paths, header guards, and a build system were added. See
["What changed" below](#what-changed-vs-the-original-files) for the exact
list of edits.

## Layout

```
ed_lanczos/
├── CMakeLists.txt          # build configuration
├── include/ed_lanczos/     # public headers (installed / included as <ed_lanczos/xxx.h>)
│   ├── common.h            # typedefs, constants, small shared utilities
│   ├── applyHamiltonian.h  # real-space Hamiltonian application (H|psi>)
│   ├── kspace_1o2.h        # spin-1/2 momentum-space (symmetry-reduced) Hamiltonian
│   ├── Lanczos.h           # real-space Lanczos ground-state solvers
│   ├── Lanczos_kspace.h    # k-space Lanczos solvers (uses ARPACK via ezarpack)
│   └── unit_cells.h        # lattice/unit-cell geometry definitions
├── src/                    # implementation (.cpp) files, one per header
├── apps/
│   └── ed_kspace_1o2_main.cpp   # `main()` driver: sets up a model and runs the solver
└── tests/                  # placeholder for a future unit-test suite
```

The convention is: `include/<project>/...` for public headers, `src/` for
implementation, `apps/` for executables with a `main()`. This is the layout
CMake's `target_include_directories` and most package managers (vcpkg,
Conan, `find_package`) expect, and it keeps "library that could be reused"
cleanly separated from "one particular study/experiment that uses it."

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
path, point CMake at them:

```bash
cmake -B build \
  -DEIGEN3_INCLUDE_DIR=/path/to/eigen \
  -DEIGENRAND_INCLUDE_DIR=/path/to/EigenRand \
  -DEZARPACK_INCLUDE_DIR=/path/to/ezarpack/include \
  -DARPACK_LIBRARY=/path/to/libarpack.so
cmake --build build -j
```

If Eigen is installed system-wide (e.g. via `apt install libeigen3-dev` or a
Conan/vcpkg toolchain file), `find_package(Eigen3)` will pick it up
automatically and you can drop `-DEIGEN3_INCLUDE_DIR`.

> **Note:** this environment does not have Eigen/EigenRand/ezarpack/ARPACK
> installed, so the reorganized project **has not been compiled or tested
> here**. The file moves and include-path edits were applied mechanically
> and reviewed by hand, but please do a test build before relying on it.

Running the driver:

```bash
./build/ed_kspace_1o2
```

## What changed vs. the original files

Purely structural/organizational changes — nothing in the numerical logic
was touched:

1. **Directory layout**: files were split into `include/ed_lanczos/`, `src/`,
   and `apps/` instead of sitting flat in one folder.
2. **Include paths**: internal includes now read `#include "ed_lanczos/common.h"`
   etc. instead of `#include "common.h"`, matching the new layout.
3. **Header guards**: `#pragma once` was added to `applyHamiltonian.h`,
   `kspace_1o2.h`, `Lanczos.h`, `Lanczos_kspace.h`, and `unit_cells.h` (only
   `common.h` had one before).
4. **`main()` isolated**: `ED_kspace_1o2.cpp` (the only file with a `main`)
   moved to `apps/ed_kspace_1o2_main.cpp`. Its local `find_index` helper was
   wrapped in an anonymous namespace so it has internal linkage instead of
   leaking an exported symbol named `find_index`.
5. **Line endings**: normalized CRLF → LF throughout.
6. **Build system**: added `CMakeLists.txt` (library target `ed_lanczos` +
   executable `ed_kspace_1o2`), `.gitignore`, and this `README.md`. None of
   these existed before.

## Recommendations not applied (would need testing against real data)

These are worth doing but are riskier changes I didn't make blind, since
there's no way to compile/test in this environment:

- `common.h` has `using namespace std;` at namespace scope in a header,
  which forces that pollution onto every translation unit that includes it
  (directly or transitively — i.e. almost the whole project). The safer
  pattern is either fully-qualifying `std::` names or scoping the directive
  to `.cpp` files only. Removing it here would touch essentially every file,
  so it's left as-is pending a real build to verify nothing breaks.
- Several `.cpp` files (`applyHamiltonian.cpp`, `kspace_1o2.cpp`,
  `unit_cells.cpp`) have little to no function-level documentation compared
  to `Lanczos_kspace.cpp`. Worth a documentation pass, especially for the
  physics conventions (bond/sublattice indexing, k-point conventions).
- `apps/ed_kspace_1o2_main.cpp` hardcodes all physical parameters (`J`, `K`,
  `Γ`, field strength, angles, lattice size, etc.) at the top of `main`.
  Moving these to a config file or command-line arguments (`argc`/`argv` are
  already accepted but unused) would make parameter sweeps much easier
  without recompiling.
- Numerous `cout <<` debug prints remain scattered through the solvers
  (`Lanczos.cpp`, `Lanczos_kspace.cpp`, `kspace_1o2.cpp`, the app). Consider
  gating these behind a verbosity flag or a small logging utility.
