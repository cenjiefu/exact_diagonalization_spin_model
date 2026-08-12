#include "ed_lanczos/Lanczos_kspace.h"

tuple<double,double,int> basicLanczos_kspace(uint32_t N_Lanczos_max, double E_tol, Eigen::Ref<SpMatCD> Hamiltonian){

	srand((unsigned int) 0);

    DVec An = DVec::Zero(N_Lanczos_max);
    DVec Bn = DVec::Zero(N_Lanczos_max);

    uint32_t phi_length = Hamiltonian.rows();
    CVec phi0 = CVec::Random(phi_length);
    phi0.normalize();
    // CVec phi_ini = phi0; // remember the initial state to calculate the eigenvectors in 2nd round

    ////
    // First round of Lanzcos for energies
    CVec phi1 = CVec::Zero(phi_length);
    phi1 = Hamiltonian*phi0;
    An(0) = (phi0.dot(phi1)).real();
    // An(0) = complex<double>(phi0.adjoint()*phi1).real();
    phi1 = phi1 - An(0)*phi0;
    Bn(0) = phi1.norm();
    phi1.normalize();

    double E1_last = 0; // energy of the ground state in the last Lanczos iter
    double E2_last = 0; // energy of the 1st excited state in the last Lanczos iter
    DVec Energy;
    DMat EV_Lanzcos;
    uint32_t N_Lanczos_used = 0;
    CVec phi2 = CVec::Zero(phi_length);
    for(uint32_t j =1; j<N_Lanczos_max; j++){
        phi2 = Hamiltonian*phi1;
        phi2 = phi2 - Bn(j-1)*phi0;
        An(j) = (phi1.dot(phi2)).real();
        // An(j) = complex<double>(phi1.adjoint()*phi2).real();
        phi2 = phi2 - An(j)*phi1;
        Bn(j) = phi2.norm();
        
        phi0 = phi1;
        phi2.normalize();
        phi1 = phi2;
        
        // calculate eigenvalues and eigenvectors after 10 iters
        if (j>10){
            N_Lanczos_used = j+1;
            Eigen::SelfAdjointEigenSolver<DMat> es(j+1);
            es.computeFromTridiagonal(An.head(j+1),Bn.head(j));
            Energy = es.eigenvalues();
            if (abs(Energy(0)-E1_last)<E_tol && abs(Energy(1)-E2_last)<E_tol){
                EV_Lanzcos = es.eigenvectors(); // Lanzcos eigenstates
                break;
            }
            E1_last = Energy(0);
            E2_last = Energy(1);
        }
    }
    
    // cout << "N_Lanczos_used = " << N_Lanczos_used << endl;
    // cout << "# state = " << Energy.rows() << endl;
    // cout << fixed << setprecision(6);
    cout << Energy.head(10) << endl;

    return {Energy(0),Energy(1), N_Lanczos_used};
}

tuple<DVec,CMat,int> basicLanczos_kspace_withState(uint32_t N_Lanczos_max, double E_tol, uint32_t N_state, Eigen::Ref<SpMatCD> Hamiltonian){

    srand((unsigned int) 0);

    DVec An = DVec::Zero(N_Lanczos_max);
    DVec Bn = DVec::Zero(N_Lanczos_max);

    uint32_t phi_length = Hamiltonian.rows();
    CVec phi0 = CVec::Random(phi_length);
    phi0.normalize();
    CVec phi_ini = CVec(phi0); // remember the initial state to calculate the eigenvectors in 2nd round

    ////
    // First round of Lanzcos for energies
    CVec phi1 = CVec::Zero(phi_length);
    phi1 = Hamiltonian*phi0;
    An(0) = (phi0.dot(phi1)).real();
    // An(0) = complex<double>(phi0.adjoint()*phi1).real();
    phi1 = phi1 - An(0)*phi0;
    Bn(0) = phi1.norm();
    phi1.normalize();


    DVec E_last = DVec::Zero(N_state); // energy of the state in the last Lanczos iter
    DVec Energy;
    DMat EV_Lanzcos;
    uint32_t N_Lanczos_used = 0;
    CVec phi2 = CVec::Zero(phi_length);
    for(uint32_t j =1; j<N_Lanczos_max; j++){
        phi2 = Hamiltonian*phi1;
        phi2 = phi2 - Bn(j-1)*phi0;
        An(j) = (phi1.dot(phi2)).real();
        // An(j) = complex<double>(phi1.adjoint()*phi2).real();
        phi2 = phi2 - An(j)*phi1;
        Bn(j) = phi2.norm();
        
        phi0 = phi1;
        phi2.normalize();
        phi1 = phi2;
        
        // calculate eigenvalues and eigenvectors after 10 iters
        if (j>10+N_state){
            N_Lanczos_used = j+1;
            Eigen::SelfAdjointEigenSolver<DMat> es(j+1);
            es.computeFromTridiagonal(An.head(j+1),Bn.head(j));
            Energy = es.eigenvalues();
            EV_Lanzcos = es.eigenvectors(); // Lanzcos eigenstates
            if ( ((Energy.head(N_state)-E_last).array().abs() < E_tol).all()){
                break;
            }
            E_last = Energy.head(N_state);
        }
    }
    
    // cout << "N_Lanczos_used = " << N_Lanczos_used << endl;
    // cout << "# state = " << Energy.rows() << endl;
    // cout << fixed << setprecision(6);
    // cout << Energy.head(10) << endl;
    CMat states = CMat::Zero(phi_length,N_state);
    //// Second round of Lanzcos for eigenstates
    phi0 = phi_ini;
    phi1.setZero(phi_length);
    phi1 = Hamiltonian*phi0;
    An(0) = (phi0.dot(phi1)).real();
    // An(0) = complex<double>(phi0.adjoint()*phi1).real();
    phi1 = phi1 - An(0)*phi0;
    Bn(0) = phi1.norm();
    phi1.normalize();
    states = states + phi0*EV_Lanzcos.row(0).head(N_state);

    for(uint32_t j =1; j<N_Lanczos_used; j++){
        phi2 = Hamiltonian*phi1;
        phi2 = phi2 - Bn(j-1)*phi0;
        An(j) = (phi1.dot(phi2)).real();
        // An(j) = complex<double>(phi1.adjoint()*phi2).real();
        phi2 = phi2 - An(j)*phi1;
        Bn(j) = phi2.norm();
        
        phi0 = phi1;
        phi2.normalize();
        phi1 = phi2;       
        states = states + phi0*EV_Lanzcos.row(j).head(N_state);
    }

    return {Energy.head(N_state),states, N_Lanczos_used};
}


tuple<DVec,CMat,int> fullLanczos_kspace_withState(uint32_t N_Lanczos_max, double E_tol, uint32_t N_state, Eigen::Ref<SpMatCD> Hamiltonian){

    srand((unsigned int) 0);

    DVec An = DVec::Zero(N_Lanczos_max);
    DVec Bn = DVec::Zero(N_Lanczos_max);

    uint32_t phi_length = Hamiltonian.rows();
    CMat ortho_phi = CMat::Zero(phi_length,N_Lanczos_max+1);

    ortho_phi.col(0) = CVec::Random(phi_length);
    // phi0.array() += 1;
    // phi0 /= 2;
    ortho_phi.col(0).normalize();
    
    ortho_phi.col(1) = CVec::Zero(phi_length);
    ortho_phi.col(1) = Hamiltonian*ortho_phi.col(0);
    An(0) = (ortho_phi.col(0).dot(ortho_phi.col(1))).real();
    ortho_phi.col(1) = ortho_phi.col(1) - An(0)*ortho_phi.col(0);
    Bn(0) = ortho_phi.col(1).norm();
    ortho_phi.col(1).normalize();
    // cout << "here1" << endl;

    DVec E_last = DVec::Zero(N_state); // energy of the state in the last Lanczos iter
    DVec Energy;
    DMat EV_Lanzcos;
    uint32_t N_Lanczos_used = 0;
    CVec proj;
    for(uint32_t j =1; j<N_Lanczos_max; j++){
        ortho_phi.col(j+1) = Hamiltonian*ortho_phi.col(j) - Bn(j-1)*ortho_phi.col(j-1);
        An(j) = (ortho_phi.col(j).dot(ortho_phi.col(j+1))).real();
        ortho_phi.col(j+1) = ortho_phi.col(j+1) - An(j)*ortho_phi.col(j);
        proj = (ortho_phi.leftCols(j+1)).adjoint()*ortho_phi.col(j+1);
        ortho_phi.col(j+1) = ortho_phi.col(j+1) - (ortho_phi.leftCols(j+1))*proj;
        Bn(j) = ortho_phi.col(j+1).norm();

        ortho_phi.col(j+1).normalize();
        // cout << "j+1 = " << j+1 << " " << ortho_phi.leftCols(j+1).cols() << endl;
        
        // calculate eigenvalues and eigenvectors after 10 iters
        if (j>10+N_state){
            N_Lanczos_used = j+1;
            Eigen::SelfAdjointEigenSolver<DMat> es(j+1);
            es.computeFromTridiagonal(An.head(j+1),Bn.head(j));
            Energy = es.eigenvalues();
            // cout << "here3" << endl;
            EV_Lanzcos = es.eigenvectors(); // Lanzcos eigenstates
            if ( ((Energy.head(N_state)-E_last).array().abs() < E_tol).all()){
                break;
            }
            E_last = Energy.head(N_state);
        }
        // cout << "N_Lanczos_used = " << N_Lanczos_used << endl;
    }
    
    // cout << "N_Lanczos_used = " << N_Lanczos_used << endl;
    // cout << "# state = " << Energy.rows() << endl;
    // cout << fixed << setprecision(6);
    // cout << Energy.head(10) << endl;
    // for(int j = 0; j<1; j++){
    //     for(int k = j+1; k<N_Lanczos_max+1; k++){
    //         cout << abs(ortho_phi.col(j).dot(ortho_phi.col(k))) << "    ";
    //     }
    //     cout << endl;
    // } 

    // CVec ground = ortho_phi.leftCols(N_Lanczos_used)*EV_Lanzcos.col(0);
    // CVec excited1 = ortho_phi.leftCols(N_Lanczos_used)*EV_Lanzcos.col(1);
    // CVec excited2 = ortho_phi.leftCols(N_Lanczos_used)*EV_Lanzcos.col(2);

    CMat states = ortho_phi.leftCols(N_Lanczos_used)*EV_Lanzcos.leftCols(N_state);

    return {Energy.head(N_state),states, N_Lanczos_used};
}


/*  calculate overlap of each Lanczos vector with initial random vector for finite temperature sampling */
tuple<DVec,CVec,DMat,DMat> fullLanczos_kspace_r0(uint32_t N_subspace, double E_tol, Eigen::Ref<SpMatCD> Hamiltonian){

    Eigen::Rand::P8_mt19937_64 urng{time(0)};

    DVec An = DVec::Zero(N_subspace);
    DVec Bn = DVec::Zero(N_subspace);

    uint32_t phi_length = Hamiltonian.rows();
    CMat ortho_phi = CMat::Zero(phi_length,N_subspace+1);

    // DMat phi0_real = Eigen::Rand::normal<DMat>(phi_length, 1, urng);
    // DMat phi0_imag = Eigen::Rand::normal<DMat>(phi_length, 1, urng);
    // ortho_phi.col(0) = phi0_real + phi0_imag*I;
    ortho_phi.col(0) = CVec::Random(phi_length);
    ortho_phi.col(0).normalize();
    
    ortho_phi.col(1) = CVec::Zero(phi_length);
    ortho_phi.col(1) = Hamiltonian*ortho_phi.col(0);
    An(0) = (ortho_phi.col(0).dot(ortho_phi.col(1))).real();
    ortho_phi.col(1) = ortho_phi.col(1) - An(0)*ortho_phi.col(0);
    Bn(0) = ortho_phi.col(1).norm();
    ortho_phi.col(1).normalize();

    DVec E_last = DVec::Zero(N_subspace); // energy of the state in the last Lanczos iter
    DVec Energy;
    DMat EV_Lanzcos;
    uint32_t N_Lanczos_used = 0;
    CVec proj;
    for(uint32_t j =1; j<N_subspace; j++){
        ortho_phi.col(j+1) = Hamiltonian*ortho_phi.col(j) - Bn(j-1)*ortho_phi.col(j-1);
        An(j) = (ortho_phi.col(j).dot(ortho_phi.col(j+1))).real();
        // cout << "An = " << (ortho_phi.col(j).dot(ortho_phi.col(j+1))) << endl;
        ortho_phi.col(j+1) = ortho_phi.col(j+1) - An(j)*ortho_phi.col(j);
        proj = (ortho_phi.leftCols(j+1)).adjoint()*ortho_phi.col(j+1);
        ortho_phi.col(j+1) = ortho_phi.col(j+1) - (ortho_phi.leftCols(j+1))*proj;
        Bn(j) = ortho_phi.col(j+1).norm();

        ortho_phi.col(j+1).normalize();
        
        N_Lanczos_used = j+1;
        // cout << "N_Lanczos_used = " << N_Lanczos_used << endl;
    }
    
    // calculate eigenvalues and eigenvectors
    Eigen::SelfAdjointEigenSolver<DMat> es(N_subspace);
    es.computeFromTridiagonal(An.head(N_subspace),Bn.head(N_subspace-1));
    Energy = es.eigenvalues();
    EV_Lanzcos = es.eigenvectors(); // Lanzcos eigenstates
    // ((Energy.head(N_state)-E_last).array().abs() < E_tol).all()
    // E_last = Energy.head(N_state);

    // CMat states = ortho_phi.leftCols(N_subspace)*EV_Lanzcos.leftCols(N_subspace);
    CVec overlap = EV_Lanzcos.row(0);

    DMat H_Lanzcos = DMat::Zero(N_subspace,N_subspace);
    for(uint32_t j =0; j<N_subspace-1; j++){
        H_Lanzcos(j,j) = An(j);
        H_Lanzcos(j+1,j) = Bn(j);
        H_Lanzcos(j,j+1) = Bn(j);
    }
    H_Lanzcos(N_subspace-1,N_subspace-1) = An(N_subspace-1);

    return {Energy, overlap, H_Lanzcos, EV_Lanzcos};
}

/*  calculate overlap of each Lanczos vector with initial random vector for finite temperature sampling */
tuple<DMat,CMat> fullLanczos_kspace_r0(uint32_t N_subspace, double E_tol, Eigen::Ref<SpMatCD> Hamiltonian, Eigen::Ref<CVec> r0){


    DVec An = DVec::Zero(N_subspace);
    DVec Bn = DVec::Zero(N_subspace);

    uint32_t phi_length = Hamiltonian.rows();
    CMat ortho_phi = CMat::Zero(phi_length,N_subspace+1);
    ortho_phi.col(0) = r0;
    ortho_phi.col(0).normalize();
    
    ortho_phi.col(1) = CVec::Zero(phi_length);
    ortho_phi.col(1) = Hamiltonian*ortho_phi.col(0);
    An(0) = (ortho_phi.col(0).dot(ortho_phi.col(1))).real();
    ortho_phi.col(1) = ortho_phi.col(1) - An(0)*ortho_phi.col(0);
    Bn(0) = ortho_phi.col(1).norm();
    ortho_phi.col(1).normalize();

    DVec E_last = DVec::Zero(N_subspace); // energy of the state in the last Lanczos iter
    DVec Energy;
    DMat EV_Lanzcos;
    uint32_t N_Lanczos_used = 0;
    CVec proj;
    for(uint32_t j =1; j<N_subspace; j++){
        ortho_phi.col(j+1) = Hamiltonian*ortho_phi.col(j) - Bn(j-1)*ortho_phi.col(j-1);
        An(j) = (ortho_phi.col(j).dot(ortho_phi.col(j+1))).real();
        // cout << "An = " << (ortho_phi.col(j).dot(ortho_phi.col(j+1))) << endl;
        ortho_phi.col(j+1) = ortho_phi.col(j+1) - An(j)*ortho_phi.col(j);
        proj = (ortho_phi.leftCols(j+1)).adjoint()*ortho_phi.col(j+1);
        ortho_phi.col(j+1) = ortho_phi.col(j+1) - (ortho_phi.leftCols(j+1))*proj;
        Bn(j) = ortho_phi.col(j+1).norm();

        ortho_phi.col(j+1).normalize();
        
        N_Lanczos_used = j+1;
        // cout << "N_Lanczos_used = " << N_Lanczos_used << endl;
    }
    
    // calculate eigenvalues and eigenvectors
    // Eigen::SelfAdjointEigenSolver<DMat> es(N_subspace);
    // es.computeFromTridiagonal(An.head(N_subspace),Bn.head(N_subspace-1));
    // Energy = es.eigenvalues();
    // ((Energy.head(N_state)-E_last).array().abs() < E_tol).all()
    // E_last = Energy.head(N_state);


    DMat H_Lanzcos = DMat::Zero(N_subspace+1,N_subspace+1);
    for(uint32_t j =0; j<N_subspace; j++){
        H_Lanzcos(j,j) = An(j);
        H_Lanzcos(j+1,j) = Bn(j);
        H_Lanzcos(j,j+1) = Bn(j);
    }

    return {H_Lanzcos, ortho_phi};
}
