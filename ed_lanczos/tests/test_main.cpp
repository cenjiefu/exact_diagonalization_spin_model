#include "ed_lanczos/common.h"
#include "ed_lanczos/Lanczos_kspace.h"
#include "ed_lanczos/kspace_1o2.h"
#include "ed_lanczos/applyHamiltonian.h"
#include "ed_lanczos/unit_cells.h"
#include <gtest/gtest.h>

TEST(SampleTest, Groundstate_check) {
    // Hamiltonian parameter
    double J = -1.0;
    double K = -0.5;
    double G = 0.5;
    double Gp = 0.0; 
    double Ac = 0.0;
    double Dc = 0.0;
    double J3 = 0.0;
    double ga = 1.0;
    double gc = 1.0;
    double h = 2.0;
    double th = pi/2;
    double phi = 0;

    double Jpp = K/3.0 + 2.0*(G-Gp)/3.0;
    double Jcp = K/3.0 - (G-Gp)/3.0;
    double Jab = J - Gp +Jcp;
    double Jc = J + 2.0*Gp + Jpp;


    Mat3d Dx, Dy, Dz, Dm;
    Dz << Jab+Jpp, 0, -sqrt(2)*Jcp, 
        0, Jab-Jpp, 0,
        -sqrt(2)*Jcp,    0, Jc;

    Dx << Jab-Jpp/2.0, -sqrt(3)/2.0*Jpp, sqrt(2)/2.0*Jcp,
        -sqrt(3)/2.0*Jpp, Jab+Jpp/2.0, -sqrt(6)/2.0*Jcp,
        sqrt(2)/2.0*Jcp, -sqrt(6)/2.0*Jcp, Jc;

    Dy << Jab-Jpp/2.0, sqrt(3)/2.0*Jpp, sqrt(2)/2.0*Jcp,
        sqrt(3)/2.0*Jpp, Jab+Jpp/2.0, sqrt(6)/2.0*Jcp,
        sqrt(2)/2.0*Jcp, sqrt(6)/2.0*Jcp, Jc;


    Dm << 0, Dc, 0,
        -Dc, 0, 0,
        0, 0, 0;

    string cluster_type = "site6";
    if (cluster_type=="site6")
        Dm = 3*Dm;

    auto [N_site, x_bonds, y_bonds, z_bonds, DM_bonds, J3_bonds, nx, ny, Tx_bonds, Ty_bonds] = get_unit_cell_params(cluster_type);
    string filename = cluster_type+"_1o2_h_"+to_string(h)+"_th_"+to_string(th)+"_phi_"+to_string(phi)+"_J_"+to_string(J)+"_K_"+to_string(K)+"_G_"+to_string(G)
            +"_Gp_"+to_string(Gp)+"_J3_"+to_string(J3);

    

    double ha = ga*h*sin(th)*cos(phi);
    double hb = ga*h*sin(th)*sin(phi);
    double hc = gc*h*cos(th);

    // Bundle the couplings and bond geometry
    HamiltonianParams model;
    model.N_site = N_site;
    model.Dx = Dx; model.Dy = Dy; model.Dz = Dz; model.Dm = Dm;
    model.ha = ha; model.hb = hb; model.hc = hc; model.Ac = Ac;
    model.J3 = J3;
    model.x_bonds = x_bonds; model.y_bonds = y_bonds; model.z_bonds = z_bonds;
    model.DM_bonds = DM_bonds; model.J3_bonds = J3_bonds;

    uint32_t N_Lanczos_max = 500;
    double E_tol = 1e-8; // energy convergence tolerance
    uint32_t N_state = 4;

    // Print ezARPACK version
    std::cout << "Using ezARPACK version " << EZARPACK_VERSION << std::endl;
    using solver_t = ezarpack::arpack_solver<ezarpack::Complex, ezarpack::eigen_storage>;
    bool EV_flag = false;
    
    // Specify parameters for the solver
    using params_t = solver_t::params_t;
    ////     Compute Ritz vectors (eigenvectors).
    params_t::compute_vectors_t compute_vectors;
    if (EV_flag)
        compute_vectors = params_t::Ritz;
    else
        compute_vectors = params_t::None;
        
    params_t params(N_state,               // Number of low-lying eigenvalues
                  params_t::SmallestReal, // We want the smallest eigenvalues
                  compute_vectors);              // Yes, we want the eigenvectors
                                      // (Ritz vectors) as well
    params.tolerance = E_tol;
    params.max_iter = N_Lanczos_max;

    SpMatCD Hamiltonian;
    vector<vector<uint32_t>> k_tables(nx*ny);
    vector<vector<uint32_t>> k_tables_period(nx*ny);
    DMat energy1(nx*ny,N_state); // save the energy of the first 5 states for each k
    DMat energy2(nx*ny,N_state);
    DMat energy3(nx*ny,N_state); 
    vector<CMat> states(nx*ny); // kspace states for each k
    vector<CMat> states_real(nx*ny); // real states for each k
    DVec kx_list = DVec::Zero(nx*ny);
    DVec ky_list = DVec::Zero(nx*ny);
    DVec N_Lan1 = DVec::Zero(nx*ny);

    DVec energy_realspace(uint32_t(pow(2,N_site)));
    uint32_t index = 0;
    
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed_time;
    start = std::chrono::system_clock::now();

    for(uint32_t j = 0; j<nx*ny; j++){
        uint32_t kx = j%nx;
        uint32_t ky = j/nx;
        kx_list(j) = double(kx)/double(nx);
        ky_list(j) = double(ky)/double(ny);

        KSpaceGeometry geom;
        geom.N_site = N_site; geom.kx = kx; geom.ky = ky;
        geom.Tx_bonds = Tx_bonds; geom.Ty_bonds = Ty_bonds;
        geom.nx = nx; geom.ny = ny;

        cout << "# Building K table... (" << elapsed_time.count() << "s)" << endl;
        tie(k_tables[j],k_tables_period[j]) = buildKTable_1o2(geom);
        end = std::chrono::system_clock::now();
        elapsed_time = end-start;
        cout << "# Building Sparse Hamiltonian... (" << elapsed_time.count() << "s)" << endl;
        Hamiltonian = buildSparseHamiltonian_1o2(k_tables[j], k_tables_period[j], model, geom);
        end = std::chrono::system_clock::now();
        elapsed_time = end-start;
        cout << "# Solving Hamiltonian... (" << elapsed_time.count() << "s)" << endl;
        

        // First method: simple Lanczos method with full reorthogonalization
        DVec E_temp;
        tie(E_temp,states[j],N_Lan1(j)) = fullLanczos_kspace_withState(N_Lanczos_max,E_tol,N_state,Hamiltonian);
        end = std::chrono::system_clock::now();
        elapsed_time = end-start;
        energy1.row(j) = E_temp.transpose();
        cout << "kx = "+to_string(kx)+", ky = "+to_string(ky)<< endl << E_temp << endl;

        if(EV_flag){
            cout << "# Converting to real space states... (" << elapsed_time.count() << "s)" << endl;
            CMat temp = CMat::Zero(pow(2,N_site),5);
            for(uint32_t k = 0; k<N_state; k++){
               temp.col(k) = k2RealSpaceState_1o2(states[j].col(k), k_tables[j], k_tables_period[j], geom);
            }
            states_real[j] = temp;

            cout << "j = " << j << ", k space norm = " << states[j].col(0).dot(states[j].col(0)) << endl;
            cout << "j = " << j << ", real space norm = " << states_real[j].col(0).dot(states_real[j].col(0)) << endl;
        }


        // Second method: implicitly restarted Arnoldi method in the arpack-ng package
        solver_t solver(Hamiltonian.rows());
        // Linear operator representing multiplication of a given vector by our matrix
        // The operator must act on the 'in' vector and store results in 'out'.
        auto matrix_op = [Hamiltonian](solver_t::vector_const_view_t in, solver_t::vector_view_t out) {
            out = Hamiltonian*in;
        };
        // Run diagonalization!
        solver(matrix_op, params);
        end = std::chrono::system_clock::now();
        elapsed_time = end-start;
        cout << "Converged. Elapsed time in seconds: " << elapsed_time.count() << " sec" << endl;
        auto const& lambda = solver.eigenvalues();
        energy2.row(j) = lambda.real().reverse();
        // Number of converged eigenvalues
        std::cout << solver.nconv() << " out of " << params.n_eigenvalues
                << " eigenvalues have converged" << std::endl;
        // Print found eigenvalues
        std::cout << "Eigenvalues (Ritz values):" << std::endl;
        std::cout << lambda.head(N_state) << std::endl;
        // Print some computation statistics
        auto stats = solver.stats();
        std::cout << "Number of Lanczos update iterations: " << stats.n_iter << std::endl;
        N_Lan1(j) = stats.n_iter;
        std::cout << "Total number of OP*x operations: " << stats.n_op_x_operations << std::endl;
        std::cout << "Total number of steps of re-orthogonalization: " << stats.n_reorth_steps << std::endl;

        if(EV_flag){
            auto const& v = solver.eigenvectors();
            states[j] = v.colwise().reverse();
            cout << "# Converting to real space states... (" << elapsed_time.count() << "s)" << endl;
            CMat temp = CMat::Zero(pow(2,N_site),N_state);
            for(uint32_t k = 0; k<N_state; k++){
               temp.col(k) = k2RealSpaceState_1o2(states[j].col(k), k_tables[j], k_tables_period[j], geom);
            }
            states_real[j] = temp;
        }


        // Third method: full matrix diagonalization, only feasible for 6 or 12 sites 
        CMat Hamil = Eigen::MatrixXcd(Hamiltonian);
        Eigen::SelfAdjointEigenSolver<CMat> es(Hamil);
        energy3.row(j) = es.eigenvalues().head(N_state);
        cout << "kx = "+to_string(kx)+", ky = "+to_string(ky)<< endl << es.eigenvalues().head(N_state) << endl;
        energy_realspace.segment(index,k_tables[j].size()) = es.eigenvalues();
        index += k_tables[j].size();
    }


    EXPECT_TRUE((energy1-energy3).norm() < E_tol*10 and (energy2-energy3).norm() < E_tol*10);
}

int main(int argc, char *argv[])
{

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();


}

