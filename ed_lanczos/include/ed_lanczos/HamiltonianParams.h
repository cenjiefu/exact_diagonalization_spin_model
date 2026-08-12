#pragma once
#include "ed_lanczos/common.h"

// Bundles the couplings and bond geometry that define a spin Hamiltonian:
//
//   H = sum_{<ij> in x-bonds}  S_i . Dx . S_j
//     + sum_{<ij> in y-bonds}  S_i . Dy . S_j
//     + sum_{<ij> in z-bonds}  S_i . Dz . S_j
//     + sum_{<ij> in DM-bonds} DM-type coupling with tensor Dm
//     + Zeeman field (ha, hb, hc) + single-ion anisotropy Ac
//     + a further-neighbor J3 coupling on J3_bonds

struct HamiltonianParams {
    uint32_t N_site = 0;

    Mat3d Dx = Mat3d::Zero();
    Mat3d Dy = Mat3d::Zero();
    Mat3d Dz = Mat3d::Zero();
    Mat3d Dm = Mat3d::Zero();   // DM-bond exchange tensor

    double ha = 0.0, hb = 0.0, hc = 0.0; // Zeeman field components
    double Ac = 0.0;                     // single-ion anisotropy
    double J3 = 0.0;

    ArrXXi x_bonds, y_bonds, z_bonds; // bond index arrays, one row per bond
    ArrXXi DM_bonds, J3_bonds;

};
