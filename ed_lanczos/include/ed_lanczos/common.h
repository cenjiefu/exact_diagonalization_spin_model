#pragma once
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <chrono>
#include <functional>
#include <cmath>

#include <string>
#include <iterator>
#include <vector>
#include <algorithm>
#include <complex>
#include <tuple>
#include <set>

#include <omp.h>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "EigenRand/EigenRand"
#include <ezarpack/arpack_solver.hpp>
#include <ezarpack/storages/eigen.hpp>
#include <ezarpack/version.hpp>

// #include <slepceps.h>
// #include "armadillo-10.7.1/include/armadillo"

using namespace std;

const double pi = 3.14159265358979323846264338327950288419L;
const std::complex<double> I = std::complex<double>(0.0,1.0);
const double etol = 1e-8; // eigenvalue tolerance
// const double vtol = 1e-12; // eigenvector tolerance
const double vtol = 1e-7; // eigenvector tolerance
const double tol = 1e-6; // general tolerance criterion

typedef Eigen::Vector3d Vec3;
typedef Eigen::Vector4d Vec4;
typedef Eigen::Matrix3d Mat3d;
typedef Eigen::ArrayXXi ArrXXi;
// typedef Eigen::Vector<char,Dynamic> charVec;

typedef Eigen::SparseMatrix< std::complex<double>, Eigen::RowMajor > SpMatCD;
typedef Eigen::Triplet< std::complex<double> > TripletCD;
typedef Eigen::VectorXcd CVec;
typedef Eigen::VectorXd DVec;
typedef Eigen::MatrixXd DMat;
typedef Eigen::MatrixXcd CMat;
typedef Eigen::Matrix<long double, Eigen::Dynamic, 1> LDVec;
typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> LDMat;

typedef Eigen::SparseMatrix< std::complex<float>, Eigen::RowMajor > SpMatCF; 
typedef Eigen::Triplet< std::complex<float> > TripletCF;
typedef Eigen::VectorXcf CFVec;
typedef Eigen::VectorXf FVec;

typedef ezarpack::arpack_solver<ezarpack::Symmetric, ezarpack::eigen_storage>::vector_const_view_t Vec_d_const;
typedef ezarpack::arpack_solver<ezarpack::Symmetric, ezarpack::eigen_storage>::vector_view_t Vec_d;
typedef ezarpack::arpack_solver<ezarpack::Complex, ezarpack::eigen_storage>::vector_const_view_t CVec_d_const;
typedef ezarpack::arpack_solver<ezarpack::Complex, ezarpack::eigen_storage>::vector_view_t CVec_d;

/* Common Functions */

double brents_fun(std::function<double (double)> f, double lower, double upper, double tol, unsigned int max_iter);
int64_t findIndexInKTable(uint32_t representative_state, vector<uint32_t> & k_table);
uint32_t qua2dec(vector<unsigned char>& qua);
vector<unsigned char> dec2qua(uint32_t decimal, uint32_t N_digit);
void printState(vector<unsigned char>& qua);
void printState(uint32_t dec, uint32_t N_digit);
double maxAbsRowSum(Eigen::Ref<SpMatCD> A);