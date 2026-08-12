#include "ed_lanczos/common.h"
#include "ed_lanczos/Lanczos.h"
#include "ed_lanczos/applyHamiltonian.h"

tuple<double,int> basicLanczos(uint32_t N_site, uint32_t N_Lanczos_max, double E_tol, const Mat3d& Dx, const Mat3d& Dy, const Mat3d& Dz, double ha, double hb, 
	double hc, double Ac, const Mat3d& Dm, const ArrXXi& x_bonds, const ArrXXi& y_bonds, const ArrXXi& z_bonds, const ArrXXi& DM_bonds){

	srand((unsigned int) 0);

    DVec An = DVec::Zero(N_Lanczos_max);
    DVec Bn = DVec::Zero(N_Lanczos_max);

    uint32_t phi_length = pow(4,N_site);
    CVec phi0 = CVec::Random(phi_length);
    phi0.normalize();
    // CVec phi_ini = phi0; // remember the initial state to calculate the eigenvectors in 2nd round

    ////
    // First round of Lanzcos for energies
    CVec phi1 = CVec::Zero(phi_length);
    applyHamiltonian(phi0,phi1,N_site,Dx,Dy,Dz,ha,hb,hc,Ac,Dm,x_bonds,y_bonds,z_bonds,DM_bonds);
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
    for(int j =1; j<N_Lanczos_max; j++){
        phi2.setZero(phi_length);
        applyHamiltonian(phi1,phi2,N_site,Dx,Dy,Dz,ha,hb,hc,Ac,Dm,x_bonds,y_bonds,z_bonds,DM_bonds);
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
    // cout << Energy.head(10) << endl;

    return {Energy(1)-Energy(0), N_Lanczos_used};
}

tuple<double,int> basicLanczos_1o2(uint32_t N_site, uint32_t N_Lanczos_max, double E_tol, const Mat3d& Dx, const Mat3d& Dy, const Mat3d& Dz, double ha, double hb, 
    double hc, double Ac, const Mat3d& Dm, const ArrXXi& x_bonds, const ArrXXi& y_bonds, const ArrXXi& z_bonds, const ArrXXi& DM_bonds){

    srand((unsigned int) 0);

    DVec An = DVec::Zero(N_Lanczos_max);
    DVec Bn = DVec::Zero(N_Lanczos_max);

    uint32_t phi_length = pow(2,N_site);
    CVec phi0 = CVec::Random(phi_length);
    phi0.normalize();
    // CVec phi_ini = phi0; // remember the initial state to calculate the eigenvectors in 2nd round

    ////
    // First round of Lanzcos for energies
    CVec phi1 = CVec::Zero(phi_length);
    applyHamiltonian_1o2(phi0,phi1,N_site,Dx,Dy,Dz,ha,hb,hc,Ac,Dm,x_bonds,y_bonds,z_bonds,DM_bonds);
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
    for(int j =1; j<N_Lanczos_max; j++){
        phi2.setZero(phi_length);
        applyHamiltonian_1o2(phi1,phi2,N_site,Dx,Dy,Dz,ha,hb,hc,Ac,Dm,x_bonds,y_bonds,z_bonds,DM_bonds);
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
    // cout << Energy.head(10) << endl;

    return {Energy(1)-Energy(0), N_Lanczos_used};
}

tuple<double,CVec,CVec,CVec,int> basicLanczos_1o2_withState(uint32_t N_site, uint32_t N_Lanczos_max, double E_tol, const Mat3d& Dx, const Mat3d& Dy, const Mat3d& Dz, double ha, double hb, 
    double hc, double Ac, const Mat3d& Dm, const ArrXXi& x_bonds, const ArrXXi& y_bonds, const ArrXXi& z_bonds, const ArrXXi& DM_bonds){

    srand((unsigned int) 0);

    DVec An = DVec::Zero(N_Lanczos_max);
    DVec Bn = DVec::Zero(N_Lanczos_max);

    uint32_t phi_length = pow(2,N_site);
    CVec phi0 = CVec::Random(phi_length);
    // phi0.array() += 1;
    // phi0 /= 2;
    phi0.normalize();
    CVec phi_ini = phi0; // remember the initial state to calculate the eigenvectors in 2nd round
    
    ////
    // First round of Lanzcos for energies
    CVec phi1 = CVec::Zero(phi_length);
    applyHamiltonian_1o2(phi0,phi1,N_site,Dx,Dy,Dz,ha,hb,hc,Ac,Dm,x_bonds,y_bonds,z_bonds,DM_bonds);
    An(0) = (phi0.dot(phi1)).real();
    // An(0) = complex<double>(phi0.adjoint()*phi1).real();
    phi1 = phi1 - An(0)*phi0;
    Bn(0) = phi1.norm();
    phi1.normalize();

    double E0_last = 0; // energy of the ground state in the last Lanczos iter
    double E1_last = 0; // energy of the 1st excited state in the last Lanczos iter
    double E2_last = 0;
    DVec Energy;
    DMat EV_Lanzcos;
    uint32_t N_Lanczos_used = 0;
    CVec phi2 = CVec::Zero(phi_length);
    for(int j =1; j<N_Lanczos_max; j++){
        phi2.setZero(phi_length);
        applyHamiltonian_1o2(phi1,phi2,N_site,Dx,Dy,Dz,ha,hb,hc,Ac,Dm,x_bonds,y_bonds,z_bonds,DM_bonds);
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
            if (abs(Energy(0)-E0_last)<E_tol && abs(Energy(1)-E1_last)<E_tol && abs(Energy(2)-E2_last)<E_tol){
                EV_Lanzcos = es.eigenvectors(); // Lanzcos eigenstates
                break;
            }
            E0_last = Energy(0);
            E1_last = Energy(1);
            E2_last = Energy(2);
        }
    }
    
    // cout << "N_Lanczos_used = " << N_Lanczos_used << endl;
    // cout << "# state = " << Energy.rows() << endl;
    // cout << fixed << setprecision(6);
    cout << Energy.head(10) << endl;
    // cout << EV_Lanzcos.rows() << " " << EV_Lanzcos.cols() << endl;

    CVec ground = CVec::Zero(phi_length);
    CVec excited1 = CVec::Zero(phi_length);
    CVec excited2 = CVec::Zero(phi_length);
    //// Second round of Lanzcos for eigenstates
    phi0 = phi_ini;
    phi1.setZero(phi_length);
    applyHamiltonian_1o2(phi0,phi1,N_site,Dx,Dy,Dz,ha,hb,hc,Ac,Dm,x_bonds,y_bonds,z_bonds,DM_bonds);
    An(0) = (phi0.dot(phi1)).real();
    // An(0) = complex<double>(phi0.adjoint()*phi1).real();
    phi1 = phi1 - An(0)*phi0;
    Bn(0) = phi1.norm();
    phi1.normalize();
    ground = ground + EV_Lanzcos(0,0)*phi0;
    excited1 = excited1 + EV_Lanzcos(0,1)*phi0;
    excited2 = excited2 + EV_Lanzcos(0,2)*phi0;

    for(int j =1; j<N_Lanczos_used; j++){
        phi2.setZero(phi_length);
        applyHamiltonian_1o2(phi1,phi2,N_site,Dx,Dy,Dz,ha,hb,hc,Ac,Dm,x_bonds,y_bonds,z_bonds,DM_bonds);
        phi2 = phi2 - Bn(j-1)*phi0;
        An(j) = (phi1.dot(phi2)).real();
        // An(j) = complex<double>(phi1.adjoint()*phi2).real();
        phi2 = phi2 - An(j)*phi1;
        Bn(j) = phi2.norm();
        
        phi0 = phi1;
        phi2.normalize();
        phi1 = phi2;       
        ground = ground + EV_Lanzcos(j,0)*phi0;
        excited1 = excited1 + EV_Lanzcos(j,1)*phi0;
        excited2 = excited2 + EV_Lanzcos(j,2)*phi0;
    }


    return {Energy(1)-Energy(0), ground, excited1, excited2, N_Lanczos_used};
}

tuple<double,CVec,CVec,int> basicLanczos_withState(uint32_t N_site, uint32_t N_Lanczos_max, double E_tol, const Mat3d& Dx, const Mat3d& Dy, const Mat3d& Dz, double ha, double hb, 
    double hc, double Ac, const Mat3d& Dm, const ArrXXi& x_bonds, const ArrXXi& y_bonds, const ArrXXi& z_bonds, const ArrXXi& DM_bonds){

    srand((unsigned int) 0);

    DVec An = DVec::Zero(N_Lanczos_max);
    DVec Bn = DVec::Zero(N_Lanczos_max);

    uint32_t phi_length = pow(4,N_site);
    CVec phi0 = CVec::Random(phi_length);
    // phi0.array() += 1;
    // phi0 /= 2;
    // CVec phi0 = load_csv_arma<DVec>("phi0.csv");
    phi0.normalize();
    CVec phi_ini = phi0; // remember the initial state to calculate the eigenvectors in 2nd round

    ////
    // First round of Lanzcos for energies
    CVec phi1 = CVec::Zero(phi_length);
    applyHamiltonian(phi0,phi1,N_site,Dx,Dy,Dz,ha,hb,hc,Ac,Dm,x_bonds,y_bonds,z_bonds,DM_bonds);
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
    for(int j =1; j<N_Lanczos_max; j++){
        phi2.setZero(phi_length);
        applyHamiltonian(phi1,phi2,N_site,Dx,Dy,Dz,ha,hb,hc,Ac,Dm,x_bonds,y_bonds,z_bonds,DM_bonds);
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
    // cout << Energy.head(10) << endl;
    // cout << An << endl;

    CVec ground = CVec::Zero(phi_length);
    CVec excited1 = CVec::Zero(phi_length);
    //// Second round of Lanzcos for eigenstates
    phi0 = phi_ini;
    phi1.setZero(phi_length);
    applyHamiltonian(phi0,phi1,N_site,Dx,Dy,Dz,ha,hb,hc,Ac,Dm,x_bonds,y_bonds,z_bonds,DM_bonds);
    An(0) = (phi0.dot(phi1)).real();
    // An(0) = complex<double>(phi0.adjoint()*phi1).real();
    phi1 = phi1 - An(0)*phi0;
    Bn(0) = phi1.norm();
    phi1.normalize();
    ground = ground + EV_Lanzcos(0,0)*phi0;
    excited1 = excited1 + EV_Lanzcos(0,1)*phi0;

    for(int j =1; j<N_Lanczos_used; j++){
        phi2.setZero(phi_length);
        applyHamiltonian(phi1,phi2,N_site,Dx,Dy,Dz,ha,hb,hc,Ac,Dm,x_bonds,y_bonds,z_bonds,DM_bonds);
        phi2 = phi2 - Bn(j-1)*phi0;
        An(j) = (phi1.dot(phi2)).real();
        // An(j) = complex<double>(phi1.adjoint()*phi2).real();
        phi2 = phi2 - An(j)*phi1;
        Bn(j) = phi2.norm();
        
        phi0 = phi1;
        phi2.normalize();
        phi1 = phi2;       
        ground = ground + EV_Lanzcos(j,0)*phi0;
        excited1 = excited1 + EV_Lanzcos(j,1)*phi0;
    }



    return {Energy(1)-Energy(0), ground, excited1, N_Lanczos_used};
}


tuple<double,CVec,CVec,CVec,int> fullLanczos_1o2_withState(uint32_t N_site, uint32_t N_Lanczos_max, double E_tol, const Mat3d& Dx, const Mat3d& Dy, const Mat3d& Dz, double ha, double hb, 
    double hc, double Ac, const Mat3d& Dm, const ArrXXi& x_bonds, const ArrXXi& y_bonds, const ArrXXi& z_bonds, const ArrXXi& DM_bonds){

    srand((unsigned int) 0);

    uint32_t phi_length = pow(2,N_site);
    DVec An = DVec::Zero(N_Lanczos_max);
    DVec Bn = DVec::Zero(N_Lanczos_max);
    CMat ortho_phi = CMat::Zero(phi_length,N_Lanczos_max+1);

    ortho_phi.col(0) = CVec::Random(phi_length);
    // phi0.array() += 1;
    // phi0 /= 2;
    ortho_phi.col(0).normalize();
    
    ortho_phi.col(1) = CVec::Zero(phi_length);
    applyHamiltonian_1o2(ortho_phi.col(0),ortho_phi.col(1),N_site,Dx,Dy,Dz,ha,hb,hc,Ac,Dm,x_bonds,y_bonds,z_bonds,DM_bonds);
    An(0) = (ortho_phi.col(0).dot(ortho_phi.col(1))).real();
    ortho_phi.col(1) = ortho_phi.col(1) - An(0)*ortho_phi.col(0);
    Bn(0) = ortho_phi.col(1).norm();
    ortho_phi.col(1).normalize();


    double E1_last = 0; // energy of the ground state in the last Lanczos iter
    double E2_last = 0; // energy of the 1st excited state in the last Lanczos iter
    DVec Energy;
    DMat EV_Lanzcos;
    uint32_t N_Lanczos_used = 0;
    CVec proj;
    for(int j =1; j<N_Lanczos_max; j++){
        applyHamiltonian_1o2(ortho_phi.col(j),ortho_phi.col(j+1),N_site,Dx,Dy,Dz,ha,hb,hc,Ac,Dm,x_bonds,y_bonds,z_bonds,DM_bonds);
        An(j) = (ortho_phi.col(j).dot(ortho_phi.col(j+1))).real();
        proj = (ortho_phi.leftCols(j+1)).adjoint()*ortho_phi.col(j+1);
        ortho_phi.col(j+1) = ortho_phi.col(j+1) - (ortho_phi.leftCols(j+1))*proj;
        Bn(j) = ortho_phi.col(j+1).norm();

        ortho_phi.col(j+1).normalize();
        
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
    // cout << An << endl;


    CVec ground = ortho_phi.leftCols(N_Lanczos_used)*EV_Lanzcos.col(0);
    CVec excited1 = ortho_phi.leftCols(N_Lanczos_used)*EV_Lanzcos.col(1);
    CVec excited2 = ortho_phi.leftCols(N_Lanczos_used)*EV_Lanzcos.col(2);


    return {Energy(1)-Energy(0), ground, excited1, excited2, N_Lanczos_used};
}

// template <typename M>
// M load_csv_arma (const std::string & path) {
//     arma::mat X;
//     X.load(path, arma::csv_ascii);
//     return Eigen::Map<const M>(X.memptr(), X.n_rows, X.n_cols);
// }