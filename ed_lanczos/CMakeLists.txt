cmake_minimum_required(VERSION 3.16)
project(ed_lanczos LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Project-wide settings
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

option(ED_LANCZOS_ENABLE_LTO "Enable interprocedural optimization in Release builds" OFF)

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
# Eigen, EigenRand and ezarpack are header-only. If they are not installed
# where CMake/the compiler can already find them, point CMake at their
# locations, e.g.:
#   cmake -B build -DEIGEN3_INCLUDE_DIR=/path/to/eigen \
#                  -DEIGENRAND_INCLUDE_DIR=/path/to/EigenRand \
#                  -DEZARPACK_INCLUDE_DIR=/path/to/ezarpack/include
find_package(Eigen3 3.4)
find_package(OpenMP REQUIRED)

# ARPACK (via ezarpack) needs a Fortran ARPACK implementation plus BLAS/LAPACK.
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_library(ARPACK_LIBRARY NAMES arpack arpack-ng)
if(NOT ARPACK_LIBRARY)
  message(WARNING
    "ARPACK library not found automatically. Set -DARPACK_LIBRARY=/path/to/libarpack.so "
    "(or .a) when configuring, or install arpack-ng.")
endif()


# ---------------------------------------------------------------------------
# Core library: all the physics/linear-algebra code, no main()
# ---------------------------------------------------------------------------
add_library(ed_lanczos
  src/common.cpp
  src/kspace_1o2.cpp
  src/Lanczos_kspace.cpp
  src/applyHamiltonian.cpp
  src/unit_cells.cpp
)
add_library(ed_lanczos::ed_lanczos ALIAS ed_lanczos)

target_include_directories(ed_lanczos
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_include_directories(ed_lanczos PUBLIC /path/to/eigen)
target_include_directories(ed_lanczos PUBLIC /path/to/ezarpack)

target_link_libraries(ed_lanczos
  PUBLIC
    OpenMP::OpenMP_CXX
    BLAS::BLAS
    LAPACK::LAPACK
)


if(ARPACK_LIBRARY)
  target_link_libraries(ed_lanczos PUBLIC ${ARPACK_LIBRARY})
endif()

target_compile_options(ed_lanczos PRIVATE -Wall -Wextra -Wpedantic)

if(ED_LANCZOS_ENABLE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT ipo_supported)
  if(ipo_supported)
    set_target_properties(ed_lanczos PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()
endif()

# ---------------------------------------------------------------------------
# Executable(s)
# ---------------------------------------------------------------------------
add_executable(ed_kspace_1o2 apps/ed_kspace_1o2_main.cpp)
target_link_libraries(ed_kspace_1o2 PRIVATE ed_lanczos::ed_lanczos)

# ---------------------------------------------------------------------------
# Tests (optional)
# ---------------------------------------------------------------------------
option(ED_LANCZOS_BUILD_TESTS "Build unit tests" OFF)
if(ED_LANCZOS_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
