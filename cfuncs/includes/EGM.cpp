namespace EGM {

    // forward declarations
    void last_period(sol_struct*, double*, double*);
    void EGM(sol_struct*);
    void next_period_values(sol_struct*);
    void update_Q(sol_struct*);


////////////////
// 1. profile //
////////////////

void all(size_t t, par_struct *par, sol_struct *sol){
	
    size_t level = 3;
        
        logs::solve(level,"\nEGM::all\n");

    // number of points at borrowing constraint
    size_t NBC = par->Nm-par->Na;

    // parallel
    #if ONEDIM == 0 
    #pragma omp for collapse(2)
    #endif
    for(size_t i_phat  = 0;  i_phat  < par->Nphat;  i_phat++){
    for(size_t i_xihat = 0;  i_xihat < par->Nxihat; i_xihat++){

    	// a. pack
    	sol->t = t;
        sol->phat   = par->grid_phat[i_phat];
        sol->xihat  = par->grid_xihat[i_xihat];
        sol->Phat   = exp(sol->phat);

        // b. post-decision vector
        sol->a_vec = &par->grid_a[t*par->Na + 0];

        // c. last-period
        if(t == par->T-1){

            size_t i_cell    = t;        
            size_t i_states  = index::d3(i_xihat,i_phat,0,par->Nphat,par->Nm);
            
            last_period(sol, &par->m[i_cell][i_states], &par->c[i_cell][i_states]);

            continue;

        } // last period

        // d. EGM
            
            if(par->Nphat == 1){
                
                size_t i_cell    = t;        
                size_t i_states  = index::d3(i_xihat,i_phat,NBC,par->Nphat,par->Nm);

                sol->m_vec = &par->m[i_cell][i_states];        
                sol->c_vec = &par->c[i_cell][i_states];  

            }

        EGM(sol);

        // e. constraint
        if(par->Nphat == 1){

                size_t i_cell    = t; 
                size_t i_states  = index::d3(i_xihat,i_phat,0,par->Nphat,par->Nm);
            
                sol->m_vec = &par->m[i_cell][i_states];        
                sol->c_vec = &par->c[i_cell][i_states];       

            #if ONEDIM == 1 
            #pragma omp for
            #endif 
            for(size_t i_m = 0; i_m < NBC; i_m++){
                sol->m_vec[i_m] = -par->zeta[t] + (double)i_m/double(NBC)*(sol->m_vec[NBC] + par->zeta[t]);
                sol->c_vec[i_m] = sol->m_vec[i_m] + par->zeta[t];
            }

        // f. interpolate to common grid 
        } else {

            // i. interpolate to common
            size_t i_cell    = t;        
            size_t i_states  = index::d3(i_xihat,i_phat,0,par->Nphat,par->Nm);

            interpolate::evaluate_to_common(par, sol, &par->c[i_cell][i_states]);
            
            // ii. constraint
            size_t i_a = 0;
            while(isnan(sol->Q_vec[i_a])){
                i_a++;
            }
            size_t i_m = 0;
            while(i_m < par->Nm && par->grid_m[t*par->Nm + i_m] < sol->m_vec[i_a]){
                par->c[i_cell][i_states+i_m] = par->grid_m[t*par->Nm + i_m] + par->zeta[t]; 
                i_m++;
            }

            // iii. copy grid_m
            for(size_t i_m = 0; i_m < par->Nm; i_m++){
                par->m[i_cell][i_states+i_m] = par->grid_m[t*par->Nm + i_m];  
            }
        
        }

    } }

}

void last_period(sol_struct *sol, double *m_out, double *c_out){

    size_t level = 4;

        logs::solve(level,"EGM::last_period\n");    

    // a. unpack
    auto par = sol->par;

    // b. maximum m
    double m_max = par->grid_a[sol->t*par->Na + par->Na-1];

    // c. terminal consumption functions
    #if ONEDIM == 1 
    #pragma omp for
    #endif    
    for(size_t i_m = 0; i_m < par->Nm; i_m++){

        // i. cash-on-hand
        if(par->Nphat == 1){
            m_out[i_m] = 0 + (double)i_m/double(par->Nm-1)*m_max;
        } else {
            m_out[i_m] = par->grid_m[sol->t*par->Nm + i_m];
        } 

        // ii. consumption
        c_out[i_m] = m_out[i_m];
             
    }

}


////////////
// 2. EGM //
////////////

void EGM(sol_struct *sol){

    size_t level = 4;

        logs::solve(level,"EGM::EGM\n");    

    //////////////
    // 1. setup //
    //////////////

    // a. unpack
	auto par   = sol->par;
	auto t     = sol->t;
    auto a_vec = sol->a_vec;
    auto c_vec = sol->c_vec;
    auto m_vec = sol->m_vec;
    auto Q_vec = sol->Q_vec;    

    /////////////
    // 2. loop //
    /////////////
    
    size_t Nshocks = 1;

    #if ONEDIM == 1 
    #pragma omp master
    {
    #endif

    // a. initialize
	for(size_t i_a = 0; i_a < par->Na; i_a++){
	   Q_vec[i_a] = 0.0;
	}
        
    #if ONEDIM == 1 
    }
    #endif

    // b. loop over shocks
    if(t+1 < par->TR){
        Nshocks = par->Nshocks;
    }

    #if ONEDIM == 1 
    #pragma omp barrier
    #pragma omp for
    #endif
	for(size_t i_shock = 0; i_shock < Nshocks; i_shock++){

        logs::solve(level+1,"i_shock = %d\n",i_shock);
        sol->i_shock = i_shock;

        // i. next-period states
    	next_period_values(sol);

        // ii. interpolaion
    	interpolate::evaluate(par,sol);

        // iii. update
        update_Q(sol);

    }

    ////////////////
    // 4. m and c //
    ////////////////

    #if ONEDIM == 1 
    #pragma omp for
    #endif
    for(size_t i_a = 0; i_a < par->Na; i_a++){

        // a. invert Euler-equation and endogenous m
        if(isnan(Q_vec[i_a])){
            c_vec[i_a] = 0;
        } else {
            c_vec[i_a] = utilities::inv_marg_u(Q_vec[i_a],par);            
        }

        // b. endogenous cash-on-hand
        m_vec[i_a] = a_vec[i_a] + c_vec[i_a];

        logs::solve(7,"a = %g, Q = %g, c = %g, m = %g\n",
            a_vec[i_a],Q_vec[i_a],c_vec[i_a],m_vec[i_a]);

    }

}

void next_period_values(sol_struct *sol){

    size_t level = 6;

        logs::solve(level,"EGM::next_period_values\n");

    auto par = sol->par;
    auto t   = sol->t;

    /////////////
    // working //
    /////////////

    if(t+1 < par->TR){ // working
    if(par->Nphat == 1){

        // i. shocks
        double psi = par->sigma_q_psi[t]*par->iota1_x[sol->i_shock];
        double xi  = par->sigma_xi*par->iota2_x[sol->i_shock];
        double eps = 0;        
        if(mxIsInf(par->sigma_eps) == 0){
            eps = par->sigma_eps*par->iota3_x[sol->i_shock];
        }

        // ii. permanent shock
        double log_Psi_plus = par->logG[t+1] + par->mu_psi + (psi+xi)*par->K[t][0] + (psi+eps)*par->K[t][1]; 
        sol->Psi_plus = exp(log_Psi_plus);

        // iii. transitory shock
        double trans_shock = exp(par->mu_xi+par->sigma_xi*par->iota2_x[sol->i_shock]); 

        // iv. cash-on-hand
        for(size_t i_a = 0; i_a < par->Na; i_a++){
            sol->m_plus_vec[i_a] = par->R*sol->a_vec[i_a]/sol->Psi_plus + trans_shock;       
        }

    } else if(par->Nxihat == 1){

        // i. shocks
        double psi = par->sigma_q_psi[t]*par->iota1_x[sol->i_shock];
        double xi  = par->sigma_xi*par->iota2_x[sol->i_shock];
        double eps = 0;        
        if(mxIsInf(par->sigma_eps) == 0){
            eps = par->sigma_eps*par->iota3_x[sol->i_shock];
        }

        // ii. predict
        double p_predict = par->logG[t+1]+par->alpha*sol->phat + par->mu_psi;
        double y_predict = p_predict + par->mu_xi;
        double z_predict = p_predict;

        // iii. realization
        double y_plus = y_predict + psi + xi;
        double z_plus = p_predict + psi + eps;

        // iv. update
        double Delta_1 = y_plus-y_predict;
        double Delta_2 = z_plus-z_predict;

        sol->phat_plus = p_predict + Delta_1*par->K[t][0] + Delta_2*par->K[t][1];
        
        // v. cash-on-hand
        double Phat_plus = exp(sol->phat_plus);
        double Y_plus = exp(y_plus);
        sol->Psi_plus  = Phat_plus/sol->Phat;

        for(size_t i_a = 0; i_a < par->Na; i_a++){
            double M_plus = par->R*sol->a_vec[i_a]*sol->Phat + Y_plus;
            sol->m_plus_vec[i_a] = M_plus/Phat_plus;       
        }

    } else {

        double iota1_x = par->iota1_x[sol->i_shock];
        double iota2_x = par->iota2_x[sol->i_shock];
        double iota3_x = par->iota3_x[sol->i_shock];
        double iota4_x = par->iota4_x[sol->i_shock];
        double iota5_x = par->iota5_x[sol->i_shock];

        // i. predict
        double p_predict = par->logG[t+1]+par->alpha*sol->phat + par->mu_psi;
        double y_predict = p_predict + par->mu_xi + par->omega*sol->xihat;
        double z_predict = p_predict;

        // ii. realization
        double p_VsqrtD = par->VsqrtD[t][0]*iota1_x+par->VsqrtD[t][3]*iota2_x;
        double p_plus   = p_predict + par->alpha*p_VsqrtD + par->sigma_psi*iota3_x;

        double eta_VsqrtD = par->VsqrtD[t][1]*iota1_x+par->VsqrtD[t][4]*iota2_x;
        double eta_plus   = par->omega*(sol->xihat + eta_VsqrtD) + par->sigma_xi*iota4_x;

        double y_plus = p_plus + eta_plus;
        double z_plus = p_plus + par->sigma_eps*iota5_x;

        // iii. update
        double Delta_1 = y_plus-y_predict;
        double Delta_2 = z_plus-z_predict;

        sol->phat_plus  = p_predict  + Delta_1*par->K[t][0] + Delta_2*par->K[t][3];
        sol->xihat_plus = par->mu_xi + Delta_1*par->K[t][2] + Delta_2*par->K[t][5];

        // iv. cash-on-hand
        double Phat_plus = exp(sol->phat_plus);
        double Y_plus = exp(y_plus);
        sol->Psi_plus  = Phat_plus/sol->Phat;

        for(size_t i_a = 0; i_a < par->Na; i_a++){
            double M_plus = par->R*sol->a_vec[i_a]*sol->Phat + Y_plus;
            sol->m_plus_vec[i_a] = M_plus/Phat_plus;       
        }

    }


    /////////////
    // retired //
    /////////////

    } else { // retired

        // i. income
        if(par->Nphat != 1){
            sol->phat_plus = sol->phat;
        }
        if(par->Nxihat != 1){
            sol->xihat_plus = sol->xihat;
        }
        sol->Psi_plus = 1;

        // ii. cash-on-hand
        for(size_t i_a = 0; i_a < par->Na; i_a++){
            sol->m_plus_vec[i_a] = par->R*par->grid_a[sol->t*par->Na + i_a] + par->kappa;
        }

    }

}

void update_Q(sol_struct *sol){

    size_t level = 6;

        logs::solve(level,"EGM::update_Q\n");

    auto t = sol->t;
    auto par = sol->par;

    for(size_t i_a = 0; i_a < par->Na; i_a++){
		
        double weight;
        if(t+1 < par->TR){ // working
            weight = par->iota_w[sol->i_shock];
        } else {
            weight = 1.0;
        }

        double c_plus = sol->c_plus_vec[i_a];
        double Q_shocks;
		
        if(sol->c_plus_vec[i_a] <= 0){
            Q_shocks = NAN;
        } else {
            Q_shocks = par->beta*par->R*utilities::marg_u(sol->Psi_plus*c_plus,par);            
        }

        #if ONEDIM == 1 
        #pragma omp atomic
        #endif        
		sol->Q_vec[i_a] += weight*Q_shocks;
		
		logs::solve(level+1,"weight = %g, a = %g, c_plus = %g, Q_shocks = %g [%g]\n",
            weight,sol->a_vec[i_a],c_plus,Q_shocks,sol->Q_vec[i_a]);
		
    }

}

} // namespace