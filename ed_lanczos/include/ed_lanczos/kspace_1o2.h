#pragma once
#include "ed_lanczos/common.h"
#include "ed_lanczos/HamiltonianParams.h"


// Bundles the translation-symmetry data needed to build a symmetry-reduced
// (fixed momentum kx,ky) k-space Hamiltonian
// Tx_bonds and Ty_bonds define the translation
// nx and ny are the periodicity
struct KSpaceGeometry {
    uint32_t N_site, kx, ky;
    ArrXXi Tx_bonds, Ty_bonds;
    uint32_t nx, ny;

};


SpMatCD buildSparseHamiltonian_1o2(vector<uint32_t>& k_table_state, vector<uint32_t>& k_table_state_period, const HamiltonianParams& params, const KSpaceGeometry& geom);

void addBond2Hamiltonian_1o2(vector<vector<TripletCD>>& hamiltonian_matrix_elements_2d, uint32_t j, vector<uint32_t>& k_table_state, vector<uint32_t>& k_table_state_period,
	const KSpaceGeometry& geom, uint32_t state_dec, const ArrXXi& bonds, const Mat3d& D);

tuple<vector<uint32_t>,vector<uint32_t>> buildKTable_1o2(const KSpaceGeometry& geom);

CVec k2RealSpaceState_1o2(Eigen::Ref<CVec> phi_k, vector<uint32_t>& k_table_state, vector<uint32_t>& k_table_state_period,
    const KSpaceGeometry& geom);

int32_t checkState_1o2(uint32_t state_dec, const KSpaceGeometry& geom);

tuple<uint32_t,uint32_t,uint32_t> findRepresentative_1o2(uint32_t state_dec, const KSpaceGeometry& geom);

uint32_t translate_1o2(uint32_t state_dec, const ArrXXi& T_bonds, uint32_t n_times);


