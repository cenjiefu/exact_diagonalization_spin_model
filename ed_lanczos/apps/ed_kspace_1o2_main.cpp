#include "ed_lanczos/common.h"
#include "ed_lanczos/Lanczos_kspace.h"
#include "ed_lanczos/kspace_1o2.h"
#include "ed_lanczos/applyHamiltonian.h"
#include "ed_lanczos/unit_cells.h"

namespace {

// Locates the (cell, sublattice) indices of a site label within a bond
// index array. This helper is local to this translation unit.
tuple<uint32_t, uint32_t> find_index(uint32_t elem, const ArrXXi& arr);

tuple<uint32_t, uint32_t> find_index(uint32_t elem, const ArrXXi& arr){
    for(uint32_t i =0; i<arr.rows(); i++){
        for(uint32_t j =0; j<arr.cols(); j++){
            if (arr(i,j)==elem){
                return make_tuple(i,j);
            }
        }
    }
    throw invalid_argument( "coundn't find cell indices" );
    return make_tuple(-1,-1);
}

} // namespace

int main(int argc, char *argv[])
{
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
    double E_tol = 1e-6; // energy convergence tolerance
    uint32_t N_state = 4;

    // Print ezARPACK version
    std::cout << "Using ezARPACK version " << EZARPACK_VERSION << std::endl;
    using solver_t = ezarpack::arpack_solver<ezarpack::Complex, ezarpack::eigen_storage>;
    bool EV_flag = true;
    
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
    DMat energy(nx*ny,N_state); // save the energy of the first 5 states for each k
    vector<CMat> states(nx*ny); // kspace states for each k
    vector<CMat> states_real(nx*ny); // real states for each k
    DVec kx_list = DVec::Zero(nx*ny);
    DVec ky_list = DVec::Zero(nx*ny);
    DVec N_Lan1 = DVec::Zero(nx*ny);

    
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
        energy.row(j) = lambda.real().reverse();
        // Number of converged eigenvalues
        std::cout << solver.nconv() << " out of " << params.n_eigenvalues
                << " eigenvalues have converged" << std::endl;
        // Print found eigenvalues
        std::cout << "Eigenvalues (Ritz values):" << std::endl;
        std::cout << energy.row(j) << std::endl;
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

    }


    DMat parameters{{h, th, phi, J, K, G, Gp, J3, Dc}};
    cout << "h, th, phi, J, K, G, Gp, J3, Dc" << endl;
    cout << parameters << endl << endl;

    /* calculate <g(k)|S(q)|g(0)> and <e(k)|S(q)|g(0)> at k=q */ 
    vector<CMat> e_Sq_abc_g(N_state);
    uint32_t phi_length = states_real[0].rows();
    for(uint32_t s = 0; s<N_state; s++){
        CMat Sq = CMat::Zero(nx*ny,3);
        CMat Sj = CMat::Zero(N_site,3);
        for(uint32_t i = 0; i<nx*ny; i++){
            for(uint32_t j =0; j<N_site; j++){
                CVec result = DVec::Zero(phi_length);
                applySa_1o2(states_real[0].col(0),result,N_site,j); // this function is defined in real space codes
                Sj(j,0) = (states_real[i].col(s).dot(result));

                result.setZero(phi_length);
                applySb_1o2(states_real[0].col(0),result,N_site,j);
                Sj(j,1) = (states_real[i].col(s).dot(result));

                result.setZero(phi_length);
                applySc_1o2(states_real[0].col(0),result,N_site,j);
                Sj(j,2) = (states_real[i].col(s).dot(result));
            }
            if (s==0 && i==0){
                cout << "moment direction" << endl;
                cout << Sj << endl << endl;
            }
            

            // Fourier transform
            for(uint32_t cell=0; cell<x_bonds.rows(); cell++){
                uint32_t site1 = x_bonds(cell,0)-1; 
                uint32_t site2 = x_bonds(cell,1)-1;
                uint32_t cell_x = cell%nx;
                uint32_t cell_y = cell/nx;
                uint32_t qx = i%nx;
                uint32_t qy = i/nx;
                complex<double> ab_phase = exp(2*pi*I*( ((double)(qx)/(double)nx)*(1.0) + ((double)(qy)/(double)ny)*(-1.0/3.0) ));
                Sq(i,0) += (Sj(site1,0) + ab_phase*Sj(site2,0))*exp(2*pi*I*( ((double)(qx)/(double)nx)*(double)cell_x + ((double)(qy)/(double)ny)*(double)cell_y ));
                Sq(i,1) += (Sj(site1,1) + ab_phase*Sj(site2,1))*exp(2*pi*I*( ((double)(qx)/(double)nx)*(double)cell_x + ((double)(qy)/(double)ny)*(double)cell_y ));
                Sq(i,2) += (Sj(site1,2) + ab_phase*Sj(site2,2))*exp(2*pi*I*( ((double)(qx)/(double)nx)*(double)cell_x + ((double)(qy)/(double)ny)*(double)cell_y ));

            }

        }
        e_Sq_abc_g[s] = Sq;

    }

    cout << "theta = "+to_string(th)+ " ("+ to_string(th/pi)+ " pi)" << endl;
    cout << "gap(k) = \n" << energy.col(1) - energy.col(0) << endl;
    cout << "|<g(k)|S(q)|g(0)>|^2" << endl;
    cout << e_Sq_abc_g[0].array().abs2() << endl;


    DMat data1(nx*ny, 2+N_state+N_state*3+1);
    // data << kx_list, ky_list, energy.col(0), energy.col(1), g_Sq_abc_g.array().abs2(), e_Sq_abc_g.array().abs2(), N_Lan1;
    data1.col(0) = kx_list;
    data1.col(1) = ky_list;
    data1.middleCols(2, N_state) = energy;
    for(uint32_t s = 0; s<N_state; s++){
        data1.middleCols(2+N_state+s*3, 3) = e_Sq_abc_g[s].array().abs2();
    }
    data1.middleCols(2+N_state+N_state*3, 1) = N_Lan1;

    ofstream MyFile(filename+".txt");
    MyFile << parameters << endl;
    MyFile << data1 << endl;
    MyFile << endl;


    /* calculate <g|S(k)S(0)|g> for the spin */ 
    vector<CMat> g_SqSq_g(1);
    // uint32_t phi_length = states_real[0].rows();
    CMat Sq = CMat::Zero(nx*ny,3);
    CMat Sq_orb = CMat::Zero(nx*ny,3);
    complex<double> Sij[N_site][N_site][3] = {0.0};
    for(uint32_t i =0; i<N_site; i++){
        for(uint32_t j =0; j<N_site; j++){
            CVec result = DVec::Zero(phi_length);
            CVec result2 = DVec::Zero(phi_length);
            applySa_1o2(states_real[0].col(0),result,N_site,j); // this function is defined in real space codes
            applySa_1o2(result,result2,N_site,i);
            Sij[i][j][0] = (states_real[0].col(0).dot(result2));

            result.setZero(phi_length);
            result2.setZero(phi_length);
            applySb_1o2(states_real[0].col(0),result,N_site,j);
            applySb_1o2(result,result2,N_site,i);
            Sij[i][j][1]  = (states_real[0].col(0).dot(result2));

            result.setZero(phi_length);
            result2.setZero(phi_length);
            applySc_1o2(states_real[0].col(0),result,N_site,j);
            applySc_1o2(result,result2,N_site,i);
            Sij[i][j][2] = (states_real[0].col(0).dot(result2));

        }
    }


    for(uint32_t i = 0; i<nx*ny; i++){

        // Fourier transform
        for(uint32_t s1 = 0; s1<N_site; s1++){
            for(uint32_t s2 = 0; s2<N_site; s2++){
                int cell1, sub1;
                tie(cell1, sub1) =  find_index(s1+1, x_bonds);
                int cell2, sub2;
                tie(cell2, sub2) =  find_index(s2+1, x_bonds);

                int dx = cell1%nx - cell2%nx;
                int dy = cell1/nx - cell2/nx;
                uint32_t qx = i%nx;
                uint32_t qy = i/nx;

                complex<double> ab_phase;
                if (cluster_type == "site24_t2"){
                    ab_phase = exp(double(sub1-sub2)*2*pi*I*( ((double)(qx)/(double)nx)*(2.0/3.0) + ((double)(qy)/(double)ny)*(-1.0/3.0) ));
                }else{
                    ab_phase = exp(double(sub1-sub2)*2*pi*I*( ((double)(qx)/(double)nx)*(1.0) + ((double)(qy)/(double)ny)*(-1.0/3.0) ));
                }
                
                Sq(i,0) += (Sij[s1][s2][0]*ab_phase)*exp(2*pi*I*( ((double)(qx)/(double)nx)*(double)dx + ((double)(qy)/(double)ny)*(double)dy ));
                Sq(i,1) += (Sij[s1][s2][1]*ab_phase)*exp(2*pi*I*( ((double)(qx)/(double)nx)*(double)dx + ((double)(qy)/(double)ny)*(double)dy ));
                Sq(i,2) += (Sij[s1][s2][2]*ab_phase)*exp(2*pi*I*( ((double)(qx)/(double)nx)*(double)dx + ((double)(qy)/(double)ny)*(double)dy ));
                
            }
        }

    }

    g_SqSq_g[0] = Sq;


    cout << "gap(k) = \n" << energy.col(1) - energy.col(0) << endl;
    cout << "|<g(0)|SqSq|g(0)>|^2" << endl;
    cout << "spin lattice:" << endl;
    cout << g_SqSq_g[0].array().abs2() << endl;

    int L = g_SqSq_g.size();
    DMat data2(nx*ny, 2+N_state+L*3+1);
    // data << kx_list, ky_list, energy.col(0), energy.col(1), g_Sq_abc_g.array().abs2(), e_Sq_abc_g.array().abs2(), N_Lan1;
    data2.col(0) = kx_list;
    data2.col(1) = ky_list;
    data2.middleCols(2, N_state) = energy;
    for(uint32_t s = 0; s<L; s++){
        data2.middleCols(2+N_state+s*3, 3) = g_SqSq_g[s].array().abs2();
    }
    data2.middleCols(2+N_state+L*3, 1) = N_Lan1;


    MyFile << data2 << endl;
    MyFile << endl;
    MyFile.close();


    return 0;
}

