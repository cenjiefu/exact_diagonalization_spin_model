#pragma once
#include "ed_lanczos/common.h"
#include "unsupported/Eigen/MatrixFunctions"


tuple<double,double,int> basicLanczos_kspace(uint32_t N_Lanczos_max, double E_tol, Eigen::Ref<SpMatCD> Hamiltonian);

tuple<DVec,CMat,int> basicLanczos_kspace_withState(uint32_t N_Lanczos_max, double E_tol, uint32_t N_state, Eigen::Ref<SpMatCD> Hamiltonian);

tuple<DVec,CMat,int> fullLanczos_kspace_withState(uint32_t N_Lanczos_max, double E_tol, uint32_t N_state, Eigen::Ref<SpMatCD> Hamiltonian);

tuple<DVec,CVec,DMat,DMat> fullLanczos_kspace_r0(uint32_t N_subspace, double E_tol, Eigen::Ref<SpMatCD> Hamiltonian);
tuple<DMat,CMat> fullLanczos_kspace_r0(uint32_t N_subspace, double E_tol, Eigen::Ref<SpMatCD> Hamiltonian, Eigen::Ref<CVec> r0);
