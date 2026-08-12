#include "ed_lanczos/applyHamiltonian.h"

void applyHamiltonian(Eigen::Ref<CVec> phi0, Eigen::Ref<CVec> phit, const HamiltonianParams& params){

    uint32_t phi_length = phi0.rows();
    #pragma omp parallel for
    for(uint32_t j=0; j<phi_length; j++){
        vector<unsigned char> state = dec2qua(j,params.N_site);
        vector<unsigned char> from_state;

        // Zeeman terms, note the minus sign
        for(uint32_t k=0; k<params.N_site; k++){
            phit(j) = phit(j) - (double(state[k])-1.5)*params.hc*phi0(j);
            
            if (state[k]>0){
                from_state = state;
                from_state[k] = from_state[k]-1;
                phit(j) = phit(j) - sqrt(1.5*(1.5+1.0)-(double(state[k])-1.5)*(double(state[k])-1.0-1.5))*0.5*(params.ha-I*params.hb)*phi0(qua2dec(from_state));
            }
            if (state[k]<3){
                from_state = state;
                from_state[k] = from_state[k]+1;
                phit(j) = phit(j) - sqrt(1.5*(1.5+1.0)-(double(state[k])-1.5)*(double(state[k])+1.0-1.5))*0.5*(params.ha+I*params.hb)*phi0(qua2dec(from_state));
            }
            
        }
        
        // params.Ac term
        for(uint32_t k=0; k<params.N_site; k++){
            phit(j) = phit(j) + params.Ac*pow((double(state[k])-1.5),2)*phi0(j);
        }
        
        // Bond dependent interaction

        applyBond2State(phi0,phit,j,state,params.x_bonds,params.Dx);
        applyBond2State(phi0,phit,j,state,params.y_bonds,params.Dy);
        applyBond2State(phi0,phit,j,state,params.z_bonds,params.Dz);


        // // DM interaction with Dc vector
        // // since only D(1,2)=-D(2,1)=Dc are nonzero, only S+S- and S-S+ terms are nonzero
        for(int b=0; b<params.DM_bonds.rows(); b++){
            uint32_t site1 = params.DM_bonds(b,0)-1; 
            uint32_t site2 = params.DM_bonds(b,1)-1;
            // S+S- terms
            if (state[site1]>0 && state[site2]<3){
                from_state = state;
                from_state[site1] = from_state[site1]-1;
                from_state[site2] = from_state[site2]+1;
                double m = sqrt(1.5*(1.5+1.0)-(double(state[site1])-1.5)*(double(state[site1])-1.0-1.5));
                double n = sqrt(1.5*(1.5+1.0)-(double(state[site2])-1.5)*(double(state[site2])+1.0-1.5));
                phit(j) = phit(j) + m*n*I*2.0*params.Dm(0,1)/4.0*phi0(qua2dec(from_state));
            }
            // S-S+ terms
            if (state[site1]<3 && state[site2]>0){
                from_state = state;
                from_state[site1] = from_state[site1]+1;
                from_state[site2] = from_state[site2]-1;
                double m = sqrt(1.5*(1.5+1.0)-(double(state[site1])-1.5)*(double(state[site1])+1.0-1.5));
                double n = sqrt(1.5*(1.5+1.0)-(double(state[site2])-1.5)*(double(state[site2])-1.0-1.5));
                phit(j) = phit(j) - m*n*I*2.0*params.Dm(0,1)/4.0*phi0(qua2dec(from_state));
            }
        }
    }
}

void applyBond2State(Eigen::Ref<CVec> phi0, Eigen::Ref<CVec> phit, uint32_t j, vector<unsigned char>& state, const ArrXXi& bonds, const Mat3d& D){
    for(int b=0; b<bonds.rows(); b++){
        uint32_t site1 = bonds(b,0)-1; 
        uint32_t site2 = bonds(b,1)-1;
        vector<unsigned char> from_state;
        double m, n;
        // S+S+ terms
        if (state[site1]>0 && state[site2]>0){
            from_state = state;
            from_state[site1] = from_state[site1]-1;
            from_state[site2] = from_state[site2]-1;
            m = sqrt(1.5*(1.5+1.0)-(double(state[site1])-1.5)*(double(state[site1])-1.0-1.5));
            n = sqrt(1.5*(1.5+1.0)-(double(state[site2])-1.5)*(double(state[site2])-1.0-1.5));
            phit(j) = phit(j) + m*n*(D(0,0)-D(1,1)-I*2.0*D(0,1))/4.0*phi0(qua2dec(from_state));
        }
        // S-S- terms
        if (state[site1]<3 && state[site2]<3){
            from_state = state;
            from_state[site1] = from_state[site1]+1;
            from_state[site2] = from_state[site2]+1;
            m = sqrt(1.5*(1.5+1.0)-(double(state[site1])-1.5)*(double(state[site1])+1.0-1.5));
            n = sqrt(1.5*(1.5+1.0)-(double(state[site2])-1.5)*(double(state[site2])+1.0-1.5));
            phit(j) = phit(j) + m*n*(D(0,0)-D(1,1)+I*2.0*D(0,1))/4.0*phi0(qua2dec(from_state));
        }
        // S+S- terms, note the imaginary part D(0,1)-D(1,0)=0
        if (state[site1]>0 && state[site2]<3){
            from_state = state;
            from_state[site1] = from_state[site1]-1;
            from_state[site2] = from_state[site2]+1;
            m = sqrt(1.5*(1.5+1.0)-(double(state[site1])-1.5)*(double(state[site1])-1.0-1.5));
            n = sqrt(1.5*(1.5+1.0)-(double(state[site2])-1.5)*(double(state[site2])+1.0-1.5));
            phit(j) = phit(j) + m*n*(D(0,0)+D(1,1))/4.0*phi0(qua2dec(from_state));
        }
        // S-S+ terms
        if (state[site1]<3 && state[site2]>0){
            from_state = state;
            from_state[site1] = from_state[site1]+1;
            from_state[site2] = from_state[site2]-1;
            m = sqrt(1.5*(1.5+1.0)-(double(state[site1])-1.5)*(double(state[site1])+1.0-1.5));
            n = sqrt(1.5*(1.5+1.0)-(double(state[site2])-1.5)*(double(state[site2])-1.0-1.5));
            phit(j) = phit(j) + m*n*(D(0,0)+D(1,1))/4.0*phi0(qua2dec(from_state));
        }
        // S+Sc terms
        if (state[site1]>0){
            from_state = state;
            from_state[site1] = from_state[site1]-1;
            m = sqrt(1.5*(1.5+1.0)-(double(state[site1])-1.5)*(double(state[site1])-1.0-1.5));
            n = double(state[site2])-1.5;
            phit(j) = phit(j) + m*n*(D(0,2)-I*D(1,2))/2.0*phi0(qua2dec(from_state));
        }
        // S-Sc terms
        if (state[site1]<3){
            from_state = state;
            from_state[site1] = from_state[site1]+1;
            m = sqrt(1.5*(1.5+1.0)-(double(state[site1])-1.5)*(double(state[site1])+1.0-1.5));
            n = double(state[site2])-1.5;
            phit(j) = phit(j) + m*n*(D(0,2)+I*D(1,2))/2.0*phi0(qua2dec(from_state));
        }
        // ScS+ terms
        if (state[site2]>0){
            from_state = state;
            from_state[site2] = from_state[site2]-1;
            m = double(state[site1])-1.5;
            n = sqrt(1.5*(1.5+1.0)-(double(state[site2])-1.5)*(double(state[site2])-1.0-1.5));
            phit(j) = phit(j) + m*n*(D(2,0)-I*D(2,1))/2.0*phi0(qua2dec(from_state));
        }
        // ScS- terms
        if (state[site2]<3){
            from_state = state;
            from_state[site2] = from_state[site2]+1;
            m = double(state[site1])-1.5;
            n = sqrt(1.5*(1.5+1.0)-(double(state[site2])-1.5)*(double(state[site2])+1.0-1.5));
            phit(j) = phit(j) + m*n*(D(2,0)+I*D(2,1))/2.0*phi0(qua2dec(from_state));
        }
        // SzSz terms
        m = double(state[site1])-1.5;
        n = double(state[site2])-1.5;
        phit(j) = phit(j) + m*n*D(2,2)*phi0(j);
    }
}

void applyHamiltonian_1o2(Eigen::Ref<CVec> phi0, Eigen::Ref<CVec> phit, const HamiltonianParams& params){

    uint32_t phi_length = phi0.rows();
    #pragma omp parallel for
    for(uint32_t j=0; j<phi_length; j++){
        uint64_t state = j;
        uint64_t from_state;

        // Zeeman terms, note the minus sign
        for(uint32_t k=0; k<params.N_site; k++){
            phit(j) = phit(j) - (double(((state >> k) & 1ULL))-0.5)*params.hc*phi0(j);
            
            if (((state >> k) & 1ULL)){
                from_state = state ^ (1ULL << k);
                phit(j) = phit(j) - 0.5*(params.ha-I*params.hb)*phi0(from_state);
            }
            if (!((state >> k) & 1ULL)){
                from_state = state ^ (1ULL << k);
                phit(j) = phit(j) - 0.5*(params.ha+I*params.hb)*phi0(from_state);
            }
            
        }
        
        
        // Bond dependent interaction

        applyBond2State_1o2(phi0,phit,j,state,params.x_bonds,params.Dx);
        applyBond2State_1o2(phi0,phit,j,state,params.y_bonds,params.Dy);
        applyBond2State_1o2(phi0,phit,j,state,params.z_bonds,params.Dz);


        // DM interaction with Dc vector
        // since only D(1,2)=-D(2,1)=Dc are nonzero, only S+S- and S-S+ terms are nonzero
        for(int b=0; b<params.DM_bonds.rows(); b++){
            uint32_t site1 = params.DM_bonds(b,0)-1; 
            uint32_t site2 = params.DM_bonds(b,1)-1;
            // S+S- terms
            if (!((state >> site1) & 1) && ((state >> site2) & 1)){
                from_state = state ^ (1ULL << site1);
                from_state = from_state ^ (1ULL << site2);
                phit(j) = phit(j) + I*2.0*params.Dm(0,1)/4.0*phi0((from_state));
            }
            // S-S+ terms
            if (((state >> site1) & 1) && !((state >> site2) & 1)){
                from_state = state ^ (1ULL << site1);
                from_state = from_state ^ (1ULL << site2);
                phit(j) = phit(j) - I*2.0*params.Dm(0,1)/4.0*phi0((from_state));
            }
        }
    }
}

void applyBond2State_1o2(Eigen::Ref<CVec> phi0, Eigen::Ref<CVec> phit, uint32_t j, uint32_t state, const ArrXXi& bonds, const Mat3d& D){
        for(int b=0; b<bonds.rows(); b++){
            uint32_t site1 = bonds(b,0)-1; 
            uint32_t site2 = bonds(b,1)-1;
            uint32_t from_state;
            // S+S+ terms
            if (!((state >> site1) & 1ULL) && !((state >> site2) & 1ULL)){
                from_state = state ^ (1ULL << site1);
                from_state = from_state ^ (1ULL << site2);
                phit(j) = phit(j) + (D(0,0)-D(1,1)-I*2.0*D(0,1))/4.0*phi0((from_state));
            }
            // S-S- terms
            if (((state >> site1) & 1ULL) && ((state >> site2) & 1ULL)){
                from_state = state ^ (1ULL << site1);
                from_state = from_state ^ (1ULL << site2);
                phit(j) = phit(j) + (D(0,0)-D(1,1)+I*2.0*D(0,1))/4.0*phi0((from_state));
            }
            // S+S- terms, note the imaginary part D(0,1)-D(1,0)=0
            if (!((state >> site1) & 1ULL) && ((state >> site2) & 1ULL)){
                from_state = state ^ (1ULL << site1);
                from_state = from_state ^ (1ULL << site2);
                phit(j) = phit(j) + (D(0,0)+D(1,1))/4.0*phi0((from_state));
            }
            // S-S+ terms
            if (((state >> site1) & 1ULL) && !((state >> site2) & 1ULL)){
                from_state = state ^ (1ULL << site1);
                from_state = from_state ^ (1ULL << site2);
                phit(j) = phit(j) + (D(0,0)+D(1,1))/4.0*phi0((from_state));
            }
            // S+Sc terms
            if (!((state >> site1) & 1ULL)){
                from_state = state ^ (1ULL << site1);
                phit(j) = phit(j) + (double(((state >> site2) & 1ULL))-0.5)*(D(0,2)-I*D(1,2))/2.0*phi0((from_state));
            }
            // S-Sc terms
            if (((state >> site1) & 1ULL)){
                from_state = state ^ (1ULL << site1);
                phit(j) = phit(j) + (double(((state >> site2) & 1ULL))-0.5)*(D(0,2)+I*D(1,2))/2.0*phi0((from_state));
            }
            // ScS+ terms
            if (!((state >> site2) & 1ULL)){
                from_state = state ^ (1ULL << site2);
                phit(j) = phit(j) + (double(((state >> site1) & 1ULL))-0.5)*(D(2,0)-I*D(2,1))/2.0*phi0((from_state));
            }
            // ScS- terms
            if (((state >> site2) & 1ULL)){
                from_state = state ^ (1ULL << site2);
                phit(j) = phit(j) + (double(((state >> site1) & 1ULL))-0.5)*(D(2,0)+I*D(2,1))/2.0*phi0((from_state));
            }
            // SzSz terms
            phit(j) = phit(j) + (double(((state >> site1) & 1ULL))-0.5)*(double(((state >> site2) & 1ULL))-0.5)*D(2,2)*phi0(j);
        }
}

void applySc(Eigen::Ref<CVec> phi0, CVec& phit, uint32_t N_site, uint32_t site){
    uint32_t phi_length = phi0.rows();
    #pragma omp parallel for
    for(uint32_t j=0; j<phi_length; j++){
        vector<unsigned char> state = dec2qua(j,N_site);
        phit(j) = phit(j) + (double(state[site])-1.5)*phi0(j);
    }
}

void applySc_1o2(Eigen::Ref<CVec> phi0, CVec& phit, uint32_t N_site, uint32_t site){
    uint32_t phi_length = phi0.rows();
    #pragma omp parallel for
    for(uint32_t j=0; j<phi_length; j++){
        uint32_t state = j;
        phit(j) = phit(j) + (double(((state >> site) & 1ULL))-0.5)*phi0(j);
    }
}

void applySa(Eigen::Ref<CVec> phi0, CVec& phit, uint32_t N_site, uint32_t site){
    uint32_t phi_length = phi0.rows();
    #pragma omp parallel for
    for(uint32_t j=0; j<phi_length; j++){
        vector<unsigned char> state = dec2qua(j,N_site);
        vector<unsigned char> from_state;
        if (state[site]>0){
            from_state = state;
            from_state[site] = from_state[site]-1;
            // Sa = (S+ + S-)/2
            phit(j) = phit(j) + sqrt(1.5*(1.5+1.0)-(double(state[site])-1.5)*(double(state[site])-1.0-1.5))*0.5*phi0(qua2dec(from_state));
        }
        if (state[site]<3){
            from_state = state;
            from_state[site] = from_state[site]+1;
            phit(j) = phit(j) + sqrt(1.5*(1.5+1.0)-(double(state[site])-1.5)*(double(state[site])+1.0-1.5))*0.5*phi0(qua2dec(from_state));
        }
    }
}

void applySa_1o2(Eigen::Ref<CVec> phi0, CVec& phit, uint32_t N_site, uint32_t site){
    uint32_t phi_length = phi0.rows();
    #pragma omp parallel for
    for(uint32_t j=0; j<phi_length; j++){
        uint32_t state = j;
        uint32_t from_state;
        if (((state >> site) & 1ULL)){
            from_state = state ^ (1ULL << site);
            // Sa = (S+ + S-)/2
            phit(j) = phit(j) + 0.5*phi0((from_state));
        }
        if (!((state >> site) & 1ULL)){
            from_state = state ^ (1ULL << site);
            phit(j) = phit(j) + 0.5*phi0((from_state));
        }
    }
}

void applySb(Eigen::Ref<CVec> phi0, CVec& phit, uint32_t N_site, uint32_t site){
    uint32_t phi_length = phi0.rows();
    #pragma omp parallel for
    for(uint32_t j=0; j<phi_length; j++){
        vector<unsigned char> state = dec2qua(j,N_site);
        vector<unsigned char> from_state;
        if (state[site]>0){
            from_state = state;
            from_state[site] = from_state[site]-1;
            // Sb = -1i*(S+ - S-)/2
            phit(j) = phit(j) - I*sqrt(1.5*(1.5+1.0)-(double(state[site])-1.5)*(double(state[site])-1.0-1.5))*0.5*phi0(qua2dec(from_state));
        }
        if (state[site]<3){
            from_state = state;
            from_state[site] = from_state[site]+1;
            phit(j) = phit(j) + I*sqrt(1.5*(1.5+1.0)-(double(state[site])-1.5)*(double(state[site])+1.0-1.5))*0.5*phi0(qua2dec(from_state));
        }
    }
}

void applySb_1o2(Eigen::Ref<CVec> phi0, CVec& phit, uint32_t N_site, uint32_t site){
    uint32_t phi_length = phi0.rows();
    #pragma omp parallel for
    for(uint32_t j=0; j<phi_length; j++){
        uint32_t state = j;
        uint32_t from_state;
        if (((state >> site) & 1ULL)){
            from_state = state ^ (1ULL << site);
            // Sb = -1i*(S+ - S-)/2
            phit(j) = phit(j) - I*0.5*phi0((from_state));
        }
        if (!((state >> site) & 1ULL)){
            from_state = state ^ (1ULL << site);
            phit(j) = phit(j) + I*0.5*phi0((from_state));
        }
    }
}

