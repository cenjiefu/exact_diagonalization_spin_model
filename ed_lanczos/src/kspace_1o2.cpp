#include "ed_lanczos/kspace_1o2.h"

SpMatCD buildSparseHamiltonian_1o2(vector<uint32_t>& k_table_state, vector<uint32_t>& k_table_state_period, const HamiltonianParams& params, const KSpaceGeometry& geom){

	// cout << "# Building Sparse Hamiltonian... ";
	// std::chrono::time_point<std::chrono::system_clock> start, end;
	// start = std::chrono::system_clock::now();

	uint32_t k_table_size = k_table_state.size();
	SpMatCD sparse_hamiltonian = SpMatCD(k_table_size,k_table_size);
	vector<TripletCD> temp;
	vector<vector<TripletCD>> hamiltonian_matrix_elements_2d(k_table_size, temp);
	cout << "k_table_size = " << k_table_size << endl;

	#pragma omp parallel for
	for(uint32_t j=0; j < k_table_size; j++){
		uint32_t state_dec = k_table_state[j];
        uint32_t from_state_dec;
        uint32_t rep_dec, dx, dy; // the rep state and distance to the rep state
        int64_t rep_index; // index of the rep state in k_table
        double norm, m, n;
        complex<double> phase;

        // Zeeman terms, note the minus sign in the matrix element
        for(uint32_t k=0; k<geom.N_site; k++){
        	// #pragma omp critical
         //    	hamiltonian_matrix_elements.push_back(TripletCD(j,j,-(double(state[k])-1.5)*params.hc));
         	hamiltonian_matrix_elements_2d[j].push_back(TripletCD(j,j,-(double(((state_dec >> k) & 1ULL))-0.5)*params.hc));
            
            if (((state_dec >> k) & 1ULL)){
                from_state_dec = state_dec ^ (1ULL << k);
                tie(rep_dec,dx,dy) = findRepresentative_1o2(from_state_dec,geom);
                rep_index = findIndexInKTable(rep_dec, k_table_state);
                // cout << rep_index << ", ";
                if(rep_index >= 0){
                	norm = sqrt((double)k_table_state_period[rep_index]/(double)k_table_state_period[j]);
	                phase = exp(-2*pi*I*( (double)geom.kx*(double)dx/(double)geom.nx + (double)geom.ky*(double)dy/(double)geom.ny ));
	                // #pragma omp critical
	                	// hamiltonian_matrix_elements.push_back(TripletCD(j,rep_index,-m*norm*phase*0.5*(params.ha-I*params.hb)));
	                hamiltonian_matrix_elements_2d[j].push_back(TripletCD(j,rep_index, -norm*phase*0.5*(params.ha-I*params.hb)));
                } 
                
            }
            if (!((state_dec >> k) & 1ULL)){
                from_state_dec = state_dec ^ (1ULL << k);
                tie(rep_dec,dx,dy) = findRepresentative_1o2(from_state_dec,geom);
                rep_index = findIndexInKTable(rep_dec, k_table_state); // index of the rep state in k_table
                if(rep_index >= 0){
	                norm = sqrt((double)k_table_state_period[rep_index]/(double)k_table_state_period[j]);
	                phase = exp(-2*pi*I*( (double)geom.kx*(double)dx/(double)geom.nx + (double)geom.ky*(double)dy/(double)geom.ny ));
	                // #pragma omp critical
	                // 	hamiltonian_matrix_elements.push_back(TripletCD(j,rep_index,-m*norm*phase*0.5*(params.ha+I*params.hb)));
	                hamiltonian_matrix_elements_2d[j].push_back(TripletCD(j,rep_index, -norm*phase*0.5*(params.ha+I*params.hb)));
	            }
            }
        }


        // Bond dependent interaction
        addBond2Hamiltonian_1o2(hamiltonian_matrix_elements_2d,j,k_table_state,k_table_state_period,geom,state_dec,params.x_bonds,params.Dx);
        addBond2Hamiltonian_1o2(hamiltonian_matrix_elements_2d,j,k_table_state,k_table_state_period,geom,state_dec,params.y_bonds,params.Dy);
        addBond2Hamiltonian_1o2(hamiltonian_matrix_elements_2d,j,k_table_state,k_table_state_period,geom,state_dec,params.z_bonds,params.Dz);

        // // DM interaction with Dc vector
        // // since only D(1,2)=-D(2,1)=Dc are nonzero, only S+S- and S-S+ terms are nonzero
        for(uint32_t b=0; b<params.DM_bonds.rows(); b++){
            uint32_t site1 = params.DM_bonds(b,0)-1; 
            uint32_t site2 = params.DM_bonds(b,1)-1;
            // S+S- terms
            if (!((state_dec >> site1) & 1) && ((state_dec >> site2) & 1)){
                from_state_dec = state_dec ^ (1ULL << site1);
                from_state_dec = from_state_dec ^ (1ULL << site2);
	            tie(rep_dec,dx,dy) = findRepresentative_1o2(from_state_dec,geom);
	            rep_index = findIndexInKTable(rep_dec, k_table_state);
	            if(rep_index >= 0){
	            	norm = sqrt((double)k_table_state_period[rep_index]/(double)k_table_state_period[j]);
	                phase = exp(-2*pi*I*( (double)geom.kx*(double)dx/(double)geom.nx + (double)geom.ky*(double)dy/(double)geom.ny ));
		            // #pragma omp critical
	             //    	hamiltonian_matrix_elements.push_back(TripletCD(j,rep_index,m*n*norm*phase*I*2.0*params.Dm(0,1)/4.0));
		            hamiltonian_matrix_elements_2d[j].push_back(TripletCD(j,rep_index, norm*phase*I*2.0*params.Dm(0,1)/4.0));
	            }
            }
            // S-S+ terms
            if (((state_dec >> site1) & 1) && !((state_dec >> site2) & 1)){
                from_state_dec = state_dec ^ (1ULL << site1);
                from_state_dec = from_state_dec ^ (1ULL << site2);
	            tie(rep_dec,dx,dy) = findRepresentative_1o2(from_state_dec,geom);
	            rep_index = findIndexInKTable(rep_dec, k_table_state);
	            if(rep_index >= 0){
	            	norm = sqrt((double)k_table_state_period[rep_index]/(double)k_table_state_period[j]);
	                phase = exp(-2*pi*I*( (double)geom.kx*(double)dx/(double)geom.nx + (double)geom.ky*(double)dy/(double)geom.ny ));
		            // #pragma omp critical
	             //    	hamiltonian_matrix_elements.push_back(TripletCD(j,rep_index,-m*n*norm*phase*I*2.0*params.Dm(0,1)/4.0));
		            hamiltonian_matrix_elements_2d[j].push_back(TripletCD(j,rep_index, -norm*phase*I*2.0*params.Dm(0,1)/4.0));
	            }
            }
        }

        // params.J3 terms
        // three terms: SzSz, S+S-, S-S+
        for(uint32_t b=0; b<params.J3_bonds.rows(); b++){
            uint32_t site1 = params.J3_bonds(b,0)-1; 
            uint32_t site2 = params.J3_bonds(b,1)-1;
            // S+S- terms
            if (!((state_dec >> site1) & 1) && ((state_dec >> site2) & 1)){
                from_state_dec = state_dec ^ (1ULL << site1);
                from_state_dec = from_state_dec ^ (1ULL << site2);
                tie(rep_dec,dx,dy) = findRepresentative_1o2(from_state_dec,geom);
                rep_index = findIndexInKTable(rep_dec, k_table_state);
                if(rep_index >= 0){
                    norm = sqrt((double)k_table_state_period[rep_index]/(double)k_table_state_period[j]);
                    phase = exp(-2*pi*I*( (double)geom.kx*(double)dx/(double)geom.nx + (double)geom.ky*(double)dy/(double)geom.ny ));
                    // #pragma omp critical
                 //     hamiltonian_matrix_elements.push_back(TripletCD(j,rep_index,m*n*norm*phase*I*2.0*params.Dm(0,1)/4.0));
                    hamiltonian_matrix_elements_2d[j].push_back(TripletCD(j,rep_index, norm*phase*params.J3/2.0));
                }
            }
            // S-S+ terms
            if (((state_dec >> site1) & 1) && !((state_dec >> site2) & 1)){
                from_state_dec = state_dec ^ (1ULL << site1);
                from_state_dec = from_state_dec ^ (1ULL << site2);
                tie(rep_dec,dx,dy) = findRepresentative_1o2(from_state_dec,geom);
                rep_index = findIndexInKTable(rep_dec, k_table_state);
                if(rep_index >= 0){
                    norm = sqrt((double)k_table_state_period[rep_index]/(double)k_table_state_period[j]);
                    phase = exp(-2*pi*I*( (double)geom.kx*(double)dx/(double)geom.nx + (double)geom.ky*(double)dy/(double)geom.ny ));
                    // #pragma omp critical
                 //     hamiltonian_matrix_elements.push_back(TripletCD(j,rep_index,-m*n*norm*phase*I*2.0*params.Dm(0,1)/4.0));
                    hamiltonian_matrix_elements_2d[j].push_back(TripletCD(j,rep_index, norm*phase*params.J3/2.0));
                }
            }
            // SzSz terms
            hamiltonian_matrix_elements_2d[j].push_back(TripletCD(j,j,(double(((state_dec >> site1) & 1ULL))-0.5)*(double(((state_dec >> site2) & 1ULL))-0.5)*params.J3));
        }
	}

	vector<TripletCD> hamiltonian_matrix_elements((4+3*params.x_bonds.rows()*9+params.DM_bonds.rows()*2+params.J3_bonds.rows()*3)*k_table_size);
	vector<TripletCD>::iterator it = hamiltonian_matrix_elements.begin();
	for(uint32_t j=0; j < k_table_size; j++){
	  std::copy(hamiltonian_matrix_elements_2d[j].begin(), hamiltonian_matrix_elements_2d[j].end(), it);
	  it += hamiltonian_matrix_elements_2d[j].size();
	}

	sparse_hamiltonian.setFromTriplets(hamiltonian_matrix_elements.begin(), it);
  	sparse_hamiltonian.makeCompressed();
  	// end = std::chrono::system_clock::now();
  	// std::chrono::duration<double> elapsed_time = end-start;
  	// cout << "Done. (" << elapsed_time.count() << "s)" << endl;

	return sparse_hamiltonian;
}

void addBond2Hamiltonian_1o2(vector<vector<TripletCD>>& hamiltonian_matrix_elements_2d, uint32_t j, vector<uint32_t>& k_table_state, vector<uint32_t>& k_table_state_period,
	const KSpaceGeometry& geom, uint32_t state_dec, const ArrXXi& bonds, const Mat3d& D){

    for(uint32_t b=0; b<bonds.rows(); b++){
        uint32_t site1 = bonds(b,0)-1; 
        uint32_t site2 = bonds(b,1)-1;
        uint32_t from_state_dec;
        uint32_t rep_dec, dx, dy; // the rep state and distance to the rep state
        int64_t rep_index; // index of the rep state in k_table
        double norm;
        complex<double> phase;
        double m, n;
        // S+S+ terms
        if (!((state_dec >> site1) & 1ULL) && !((state_dec >> site2) & 1ULL)){
            from_state_dec = state_dec ^ (1ULL << site1);
            from_state_dec = from_state_dec ^ (1ULL << site2);
            tie(rep_dec,dx,dy) = findRepresentative_1o2(from_state_dec,geom);
            rep_index = findIndexInKTable(rep_dec, k_table_state);
            if(rep_index >= 0){
            	norm = sqrt((double)k_table_state_period[rep_index]/(double)k_table_state_period[j]);
                phase = exp(-2*pi*I*( (double)geom.kx*(double)dx/(double)geom.nx + (double)geom.ky*(double)dy/(double)geom.ny ));
	            hamiltonian_matrix_elements_2d[j].push_back(TripletCD(j,rep_index, norm*phase*(D(0,0)-D(1,1)-I*2.0*D(0,1))/4.0));
            }
            
        }
        // S-S- terms
        if (((state_dec >> site1) & 1ULL) && ((state_dec >> site2) & 1ULL)){
            from_state_dec = state_dec ^ (1ULL << site1);
            from_state_dec = from_state_dec ^ (1ULL << site2);
            tie(rep_dec,dx,dy) = findRepresentative_1o2(from_state_dec,geom);
            rep_index = findIndexInKTable(rep_dec, k_table_state);
            if(rep_index >= 0){
            	norm = sqrt((double)k_table_state_period[rep_index]/(double)k_table_state_period[j]);
                phase = exp(-2*pi*I*( (double)geom.kx*(double)dx/(double)geom.nx + (double)geom.ky*(double)dy/(double)geom.ny ));
	            // #pragma omp critical
             //    	hamiltonian_matrix_elements.push_back(TripletCD(j,rep_index,m*n*norm*phase*(D(0,0)-D(1,1)+I*2.0*D(0,1))/4.0));
	            hamiltonian_matrix_elements_2d[j].push_back(TripletCD(j,rep_index, norm*phase*(D(0,0)-D(1,1)+I*2.0*D(0,1))/4.0));
            }
        }
        // S+S- terms, note the imaginary part D(0,1)-D(1,0)=0
        if (!((state_dec >> site1) & 1ULL) && ((state_dec >> site2) & 1ULL)){
            from_state_dec = state_dec ^ (1ULL << site1);
            from_state_dec = from_state_dec ^ (1ULL << site2);
            tie(rep_dec,dx,dy) = findRepresentative_1o2(from_state_dec,geom);
            rep_index = findIndexInKTable(rep_dec, k_table_state);
            if(rep_index >= 0){
            	norm = sqrt((double)k_table_state_period[rep_index]/(double)k_table_state_period[j]);
                phase = exp(-2*pi*I*( (double)geom.kx*(double)dx/(double)geom.nx + (double)geom.ky*(double)dy/(double)geom.ny ));
	            // #pragma omp critical
             //    	hamiltonian_matrix_elements.push_back(TripletCD(j,rep_index,m*n*norm*phase*(D(0,0)+D(1,1))/4.0));
	            hamiltonian_matrix_elements_2d[j].push_back(TripletCD(j,rep_index, norm*phase*(D(0,0)+D(1,1))/4.0));
            }
        }
        // S-S+ terms
        if (((state_dec >> site1) & 1ULL) && !((state_dec >> site2) & 1ULL)){
            from_state_dec = state_dec ^ (1ULL << site1);
            from_state_dec = from_state_dec ^ (1ULL << site2);
            tie(rep_dec,dx,dy) = findRepresentative_1o2(from_state_dec,geom);
            rep_index = findIndexInKTable(rep_dec, k_table_state);
            if(rep_index >= 0){
            	norm = sqrt((double)k_table_state_period[rep_index]/(double)k_table_state_period[j]);
                phase = exp(-2*pi*I*( (double)geom.kx*(double)dx/(double)geom.nx + (double)geom.ky*(double)dy/(double)geom.ny ));
	            // #pragma omp critical
             //    	hamiltonian_matrix_elements.push_back(TripletCD(j,rep_index,m*n*norm*phase*(D(0,0)+D(1,1))/4.0));
	            hamiltonian_matrix_elements_2d[j].push_back(TripletCD(j,rep_index, norm*phase*(D(0,0)+D(1,1))/4.0));
            }
        }
        // S+Sc terms
        if (!((state_dec >> site1) & 1ULL)){
            from_state_dec = state_dec ^ (1ULL << site1);
            tie(rep_dec,dx,dy) = findRepresentative_1o2(from_state_dec,geom);
            rep_index = findIndexInKTable(rep_dec, k_table_state);
            if(rep_index >= 0){
            	norm = sqrt((double)k_table_state_period[rep_index]/(double)k_table_state_period[j]);
                phase = exp(-2*pi*I*( (double)geom.kx*(double)dx/(double)geom.nx + (double)geom.ky*(double)dy/(double)geom.ny ));
	            // #pragma omp critical
             //    	hamiltonian_matrix_elements.push_back(TripletCD(j,rep_index,m*n*norm*phase*(D(0,2)-I*D(1,2))/2.0));
	            hamiltonian_matrix_elements_2d[j].push_back(TripletCD(j,rep_index, (double(((state_dec >> site2) & 1ULL))-0.5)*norm*phase*(D(0,2)-I*D(1,2))/2.0));
            }
        }
        // S-Sc terms
        if (((state_dec >> site1) & 1ULL)){
            from_state_dec = state_dec ^ (1ULL << site1);
            tie(rep_dec,dx,dy) = findRepresentative_1o2(from_state_dec,geom);
            rep_index = findIndexInKTable(rep_dec, k_table_state);
            if(rep_index >= 0){
            	norm = sqrt((double)k_table_state_period[rep_index]/(double)k_table_state_period[j]);
                phase = exp(-2*pi*I*( (double)geom.kx*(double)dx/(double)geom.nx + (double)geom.ky*(double)dy/(double)geom.ny ));
	            // #pragma omp critical
             //    	hamiltonian_matrix_elements.push_back(TripletCD(j,rep_index,m*n*norm*phase*(D(0,2)+I*D(1,2))/2.0));
	            hamiltonian_matrix_elements_2d[j].push_back(TripletCD(j,rep_index, (double(((state_dec >> site2) & 1ULL))-0.5)*norm*phase*(D(0,2)+I*D(1,2))/2.0));
            }
        }
        // ScS+ terms
        if (!((state_dec >> site2) & 1ULL)){
            from_state_dec = state_dec ^ (1ULL << site2);
            tie(rep_dec,dx,dy) = findRepresentative_1o2(from_state_dec,geom);
            rep_index = findIndexInKTable(rep_dec, k_table_state);
            if(rep_index >= 0){
            	norm = sqrt((double)k_table_state_period[rep_index]/(double)k_table_state_period[j]);
                phase = exp(-2*pi*I*( (double)geom.kx*(double)dx/(double)geom.nx + (double)geom.ky*(double)dy/(double)geom.ny ));
	            // #pragma omp critical
             //    	hamiltonian_matrix_elements.push_back(TripletCD(j,rep_index,m*n*norm*phase*(D(2,0)-I*D(2,1))/2.0));
	            hamiltonian_matrix_elements_2d[j].push_back(TripletCD(j,rep_index, (double(((state_dec >> site1) & 1ULL))-0.5)*norm*phase*(D(2,0)-I*D(2,1))/2.0));
            }
        }
        // ScS- terms
        if (((state_dec >> site2) & 1ULL)){
            from_state_dec = state_dec ^ (1ULL << site2);
            tie(rep_dec,dx,dy) = findRepresentative_1o2(from_state_dec,geom);
            rep_index = findIndexInKTable(rep_dec, k_table_state);
            if(rep_index >= 0){
            	norm = sqrt((double)k_table_state_period[rep_index]/(double)k_table_state_period[j]);
                phase = exp(-2*pi*I*( (double)geom.kx*(double)dx/(double)geom.nx + (double)geom.ky*(double)dy/(double)geom.ny ));
	            // #pragma omp critical
             //    	hamiltonian_matrix_elements.push_back(TripletCD(j,rep_index,m*n*norm*phase*(D(2,0)+I*D(2,1))/2.0));
            	hamiltonian_matrix_elements_2d[j].push_back(TripletCD(j,rep_index, (double(((state_dec >> site1) & 1ULL))-0.5)*norm*phase*(D(2,0)+I*D(2,1))/2.0));
            }
        }
        // SzSz terms
        // #pragma omp critical
        // 	hamiltonian_matrix_elements.push_back(TripletCD(j,j,m*n*D(2,2)));
        hamiltonian_matrix_elements_2d[j].push_back(TripletCD(j,j,(double(((state_dec >> site1) & 1ULL))-0.5)*(double(((state_dec >> site2) & 1ULL))-0.5)*D(2,2)));
    }

}

tuple<vector<uint32_t>,vector<uint32_t>> buildKTable_1o2(const KSpaceGeometry& geom){
	vector<uint32_t> k_table_state, k_table_state_period;
	uint32_t phi_length = pow(2,geom.N_site);
	for (uint32_t state_dec = 0; state_dec < phi_length; state_dec++){
		int32_t normalization = checkState_1o2(state_dec, geom);
		if (normalization > 0){
      		k_table_state.push_back(state_dec);
      		k_table_state_period.push_back(normalization);
    	}
	}
	return {k_table_state,k_table_state_period};
}

/* real-space state required to calculate <2|S(q)|1>, since S(q) doesn't commute with translation operators*/
CVec k2RealSpaceState_1o2(Eigen::Ref<CVec> phi_k, vector<uint32_t>& k_table_state, vector<uint32_t>& k_table_state_period,
    const KSpaceGeometry& geom){

    CVec phi_real = CVec::Zero(pow(2,geom.N_site));
    uint32_t k_table_size = k_table_state.size();
    #pragma omp parallel for
    for(uint32_t j=0; j < k_table_size; j++){
        uint32_t state_dec = k_table_state[j];
        uint32_t translated;
        uint32_t temp;
        for(uint32_t y=0; y<geom.ny; y++){
            translated = translate_1o2(state_dec, geom.Ty_bonds, y);
            phi_real(translated) += phi_k(j)*exp(2*pi*I*( ((double)geom.ky/(double)geom.ny)*(double)y ))/sqrt(k_table_state_period[j]);
            
            for(uint32_t x=1; x<geom.nx; x++){
                temp = 0;
                for(uint32_t b=0; b<geom.Tx_bonds.rows(); b++){
                    uint32_t site1 = geom.Tx_bonds(b,0)-1; 
                    uint32_t site2 = geom.Tx_bonds(b,1)-1;
                    temp += (1ULL << site2)*((translated >> site1) & 1);
                }
                translated = temp;
                phi_real((translated)) += phi_k(j)*exp(2*pi*I*( ((double)geom.kx/(double)geom.nx)*(double)x + ((double)geom.ky/(double)geom.ny)*(double)y ))/sqrt(k_table_state_period[j]);
            }
        }
    }
    phi_real.normalize(); // should be unnecessary
    return phi_real;
}

int32_t checkState_1o2(uint32_t state_dec, const KSpaceGeometry& geom){
	set<uint32_t> trans_set; // holds all the translations
	uint32_t translated;
    uint32_t temp;
	complex<double> phase = 0; // sum of phases
	for(uint32_t y=0; y<geom.ny; y++){
		translated = translate_1o2(state_dec, geom.Ty_bonds, y);
		trans_set.insert(translated);
		if ( translated == state_dec ){
			phase += exp(-2*pi*I*( ((double)geom.ky/(double)geom.ny)*(double)y ));
		}

		// the loop below is equivalent to translate_1o2(state_dec, geom.Tx_bonds, x) with x=1..geom.nx-1; it saves time since total number of translation ~geom.nx not ~geom.nx^2
		for(uint32_t x=1; x<geom.nx; x++){
            temp = 0;
			for(uint32_t b=0; b<geom.Tx_bonds.rows(); b++){
				uint32_t site1 = geom.Tx_bonds(b,0)-1; 
	            uint32_t site2 = geom.Tx_bonds(b,1)-1;
				temp += (1ULL << site2)*((translated >> site1) & 1);
			}
			translated = temp;
			trans_set.insert((translated));
			if ( translated == state_dec ){
				phase += exp(-2*pi*I*( ((double)geom.kx/(double)geom.nx)*(double)x + ((double)geom.ky/(double)geom.ny)*(double)y ));
			}
			// printState(translated);
		}
	}
	// set<uint32_t>::iterator itr;
	// for (itr = trans_set.begin(); itr != trans_set.end(); itr++){
 //        cout << *itr<<" ";
 //    }
 //    cout << endl;
	uint32_t min = *min_element(trans_set.begin(), trans_set.end());
	if(min != state_dec){
		return -1; // not a representative state
	}
	// return the nomalization
	// cout << "size = " << trans_set.size() << endl;
	// cout << "norm = " << norm(phase) << endl;
	return (int32_t)(trans_set.size()*norm(phase));
}

tuple<uint32_t,uint32_t,uint32_t> findRepresentative_1o2(uint32_t state_dec, 
	const KSpaceGeometry& geom){

	vector<uint32_t> trans_vec, dx, dy;
	uint32_t translated;
    uint32_t temp;
	for(uint32_t y=0; y<geom.ny; y++){
		translated = translate_1o2(state_dec, geom.Ty_bonds, y);
		trans_vec.push_back((translated));
		dx.push_back(0);
		dy.push_back(y);
		
		for(uint32_t x=1; x<geom.nx; x++){
            temp = 0;
			for(uint32_t b=0; b<geom.Tx_bonds.rows(); b++){
				uint32_t site1 = geom.Tx_bonds(b,0)-1; 
	            uint32_t site2 = geom.Tx_bonds(b,1)-1;
				temp += (1ULL << site2)*((translated >> site1) & 1);
			}
			translated = temp;
			trans_vec.push_back((translated));
			dx.push_back(x);
			dy.push_back(y);
		}
	}
	uint32_t min_index = min_element(trans_vec.begin(), trans_vec.end()) - trans_vec.begin();
	return {trans_vec[min_index], dx[min_index], dy[min_index]};
}

uint32_t translate_1o2(uint32_t state_dec, const ArrXXi& T_bonds, uint32_t n_times){
	uint32_t translated = state_dec;
	for(uint32_t j=0; j<n_times; j++){
        translated = 0;
		for(uint32_t b=0; b<T_bonds.rows(); b++){
			uint32_t site1 = T_bonds(b,0)-1; 
            uint32_t site2 = T_bonds(b,1)-1;
			translated += (1ULL << site2)*((state_dec >> site1) & 1);
		}
		state_dec = translated;
	}
	return translated;
}

