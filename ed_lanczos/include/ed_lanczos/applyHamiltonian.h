#pragma once
#include "ed_lanczos/common.h"
#include "ed_lanczos/HamiltonianParams.h"

void applyHamiltonian(Eigen::Ref<CVec> phi0, Eigen::Ref<CVec> phit, const HamiltonianParams& params);

void applyHamiltonian_1o2(Eigen::Ref<CVec> phi0, Eigen::Ref<CVec> phit, const HamiltonianParams& params);

void applyBond2State(Eigen::Ref<CVec> phi0, Eigen::Ref<CVec> phit, uint32_t j, vector<unsigned char>& state, const ArrXXi& bonds, const Mat3d& D);
void applyBond2State_1o2(Eigen::Ref<CVec> phi0, Eigen::Ref<CVec> phit, uint32_t j, uint32_t state, const ArrXXi& bonds, const Mat3d& D);

void applySa(Eigen::Ref<CVec> phi0, CVec& phit, uint32_t N_site, uint32_t site);
void applySa_1o2(Eigen::Ref<CVec> phi0, CVec& phit, uint32_t N_site, uint32_t site);
void applySb(Eigen::Ref<CVec> phi0, CVec& phit, uint32_t N_site, uint32_t site);
void applySb_1o2(Eigen::Ref<CVec> phi0, CVec& phit, uint32_t N_site, uint32_t site);
void applySc(Eigen::Ref<CVec> phi0, CVec& phit, uint32_t N_site, uint32_t site);
void applySc_1o2(Eigen::Ref<CVec> phi0, CVec& phit, uint32_t N_site, uint32_t site);

