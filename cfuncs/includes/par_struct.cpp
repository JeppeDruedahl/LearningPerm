//////////////////
// 1. variables //
//////////////////

typedef struct{

    // a. grids
    size_t Na, Nm, Nphat, Nxihat;  
    double *grid_a, *grid_m, *grid_phat, *grid_xihat;

    // b. demographics
    size_t T, TR, tmin;

    // c. preferences
    double beta, rho;

    // d. income process      
    double *logG, alpha, omega;
    double kappa, logkappa;       
    size_t Nshocks;
    double mu_psi, mu_xi;
    double *sigma_q_psi, sigma_psi, sigma_xi, sigma_eps;
    double *iota1_x, *iota2_x, *iota3_x, *iota_w;
    double **K, q0;
    double *iota4_x, *iota5_x, **VsqrtD;

    // e. assets
    double R, *zeta, A0;

    // f. simulation
    size_t simN, simT, do_MPC_MPCP;
    double sigma_eta_c, sigma_eta_y;
    double *draws_q0, *draws_psi, *draws_xi, *draws_eps;  
    double *draws_eta_y, *draws_eta_c;

    // g. EGM
    double *Q_vec;

    // output
        
        // solve
        double **c, **m;

        // simulate                                
        double *sim_p, *sim_phat;
        double *sim_xi, *sim_psi, *sim_xihat, *sim_etahat;
        double *sim_A, *sim_M;
        double *sim_MPC, *sim_MPCP, *sim_logC, *sim_logY;

        // logdiffs
        double *DlogY;
        double *DlogY_lead1;
        double *DlogY_lead2;
        double *DlogY_lead3;
        double *DlogY_lead4;
        double *DlogY_lead5;
        double *DlogY_perm;
        double *DlogY_lag1;
        double *DlogY_lag2;
        double *DlogY_lag3;
        double *DlogY_lag4;
        double *DlogY_lag5;
        double *DlogC;
        double *DlogC_lead1;
        double *DlogC_lead2;
        double *DlogC_lead3;
        double *DlogC_lead4;
        double *DlogC_lead5;      

} par_struct;


namespace par {

//////////////
// 2. setup //
//////////////

void setup(par_struct *par, mxArray *plhs[], const mxArray *prhs[], size_t type){

    ///////////////
    // 1. inputs //
    ///////////////

    if(type == 1 || type == 2){

        // a. grids   
        par->Na     = (size_t) mxGetScalar(mxGetField(prhs[0],0,"Na"));
        par->Nm     = (size_t) mxGetScalar(mxGetField(prhs[0],0,"Nm"));
        par->Nphat  = (size_t) mxGetScalar(mxGetField(prhs[0],0,"Nphat"));  
        par->Nxihat = (size_t) mxGetScalar(mxGetField(prhs[0],0,"Nxihat"));             

        par->grid_a     = (double*) mxGetPr(mxGetField(prhs[0],0,"grid_a"));
        par->grid_m     = (double*) mxGetPr(mxGetField(prhs[0],0,"grid_m"));    
        par->grid_phat  = (double*) mxGetPr(mxGetField(prhs[0],0,"grid_phat"));
        par->grid_xihat = (double*) mxGetPr(mxGetField(prhs[0],0,"grid_xihat"));           

        // b. demographics
        par->T = (size_t) mxGetScalar(mxGetField(prhs[0],0,"T"));
        par->TR = (size_t) mxGetScalar(mxGetField(prhs[0],0,"TR")); 
        par->tmin = (size_t) mxGetScalar(mxGetField(prhs[0],0,"tmin"));

        // c. preferences
        par->beta  = (double) mxGetScalar(mxGetField(prhs[0],0,"beta"));
        par->rho   = (double) mxGetScalar(mxGetField(prhs[0],0,"rho"));
        
        // d. income process
        par->logG        = (double*) mxGetPr(mxGetField(prhs[0],0,"logG"));
        par->alpha       = (double) mxGetScalar(mxGetField(prhs[0],0,"alpha"));
        par->omega       = (double) mxGetScalar(mxGetField(prhs[0],0,"omega"));    
        par->kappa       = (double) mxGetScalar(mxGetField(prhs[0],0,"kappa"));
        par->logkappa    = (double) mxGetScalar(mxGetField(prhs[0],0,"logkappa"));    
        par->Nshocks     = (size_t) mxGetScalar(mxGetField(prhs[0],0,"Nshocks"));
        par->mu_psi      = (double) mxGetScalar(mxGetField(prhs[0],0,"mu_psi"));
        par->mu_xi       = (double) mxGetScalar(mxGetField(prhs[0],0,"mu_xi"));    
        par->sigma_q_psi = (double*) mxGetPr(mxGetField(prhs[0],0,"sigma_q_psi"));
        par->sigma_psi   = (double) mxGetScalar(mxGetField(prhs[0],0,"sigma_psi"));    
        par->sigma_xi    = (double) mxGetScalar(mxGetField(prhs[0],0,"sigma_xi"));
        par->sigma_eps   = (double) mxGetScalar(mxGetField(prhs[0],0,"sigma_eps"));    
        par->iota1_x     = (double*) mxGetPr(mxGetField(prhs[0],0,"iota1_x"));
        par->iota2_x     = (double*) mxGetPr(mxGetField(prhs[0],0,"iota2_x"));
        par->iota3_x     = (double*) mxGetPr(mxGetField(prhs[0],0,"iota3_x"));
        par->iota_w      = (double*) mxGetPr(mxGetField(prhs[0],0,"iota_w"));

            if(par->Nxihat != 1){
                par->iota4_x = (double*) mxGetPr(mxGetField(prhs[0],0,"iota4_x"));
                par->iota5_x = (double*) mxGetPr(mxGetField(prhs[0],0,"iota5_x"));        
            }

        // e. assets
        par->R      = (double) mxGetScalar(mxGetField(prhs[0],0,"R"));  
        par->zeta   = (double*) mxGetPr(mxGetField(prhs[0],0,"zeta"));       
        par->A0     = (double) mxGetScalar(mxGetField(prhs[0],0,"A0"));

        #if ONEDIM == 1 
        par->Q_vec = new double[par->Na];
        #endif

    }

    // h. simulate
    if(type == 2){

        par->simT        = (size_t) mxGetScalar(mxGetField(prhs[0],0,"simT"));
        par->simN        = (size_t) mxGetScalar(mxGetField(prhs[0],0,"simN"));
        par->do_MPC_MPCP = (size_t) mxGetScalar(mxGetField(prhs[0],0,"do_MPC_MPCP"));

        par->q0 = (double) mxGetScalar(mxGetField(prhs[0],0,"q0"));
        par->sigma_eta_c = (double) mxGetScalar(mxGetField(prhs[0],0,"sigma_eta_c"));
        par->sigma_eta_y = (double) mxGetScalar(mxGetField(prhs[0],0,"sigma_eta_y"));    

        par->draws_q0    = (double*) mxGetPr(mxGetField(prhs[2],0,"q0"));
        par->draws_psi   = (double*) mxGetPr(mxGetField(prhs[2],0,"psi"));
        par->draws_xi    = (double*) mxGetPr(mxGetField(prhs[2],0,"xi"));
        par->draws_eps   = (double*) mxGetPr(mxGetField(prhs[2],0,"eps"));
        par->draws_eta_y = (double*) mxGetPr(mxGetField(prhs[2],0,"eta_y"));
        par->draws_eta_c = (double*) mxGetPr(mxGetField(prhs[2],0,"eta_c"));

    }

    // i. logdiffs
    if(type == 3){    

        par->simN = (size_t) mxGetM(mxGetField(prhs[0],0,"logY"));
        par->simT = (size_t) mxGetN(mxGetField(prhs[0],0,"logY"));

        par->sim_logY = (double*) mxGetPr(mxGetField(prhs[0],0,"logY"));
        par->sim_logC = (double*) mxGetPr(mxGetField(prhs[0],0,"logC"));        

    }


        if(type == 1){
            logs::solve(1," scalars and pointers loaded\n");        
        } else if(type == 2){
            logs::simulate(1," scalars and pointers loaded\n");                    
        } else if(type == 3){
           logs::calc_logdiffs(1," scalars and pointers loaded\n");                    
       }

    ////////////////////
    // 2. cell inputs //
    ////////////////////
    
    if(type == 1 || type == 2){

        par->K = new double*[par->T];
        auto K_cell = mxGetField(prhs[0],0,"K");
        for(size_t t = 0; t < par->T;   t++){    
            par->K[t] = (double*) mxGetPr(mxGetCell(K_cell,t));
        }
        if(par->Nxihat != 1){
            par->VsqrtD = new double*[par->T];
            auto VsqrtD_cell = mxGetField(prhs[0],0,"VsqrtD");
            for(size_t t = 0; t < par->T;   t++){    
                par->VsqrtD[t] = (double*) mxGetPr(mxGetCell(VsqrtD_cell,t));                    
            }               
        }

        if(type == 1){        
            logs::solve(1," cells loaded\n");        
        } else {
            logs::simulate(1," cells loaded\n");                    
        }

    }

    ////////////////////////
    // 3. outputs - solve //
    ////////////////////////

    if(type == 1){

        // a. struct
        const char *field_names[] = {"c", "m"};
        
        size_t num_fields = sizeof(field_names)/sizeof(*field_names);
        plhs[0] = mxCreateStructMatrix(1, 1, num_fields, field_names);
        auto sol_struct = plhs[0];

        // b. cell dimensions
        size_t ndim_cell  = 1;
        auto dims_cell = new size_t[1];

        // c. array dimensions
        auto ndim = new size_t[1];
        auto dims = new size_t*[1];
                        
        // d. solution

            // cell
            dims_cell[0] = par->T;

            // array      
            ndim[0] = 3;
            dims[0] = new size_t[3];           
            dims[0][0] = par->Nm;
            dims[0][1] = par->Nphat;
            dims[0][2] = par->Nxihat;

        par->c = mymex::set_field_cell(sol_struct,"c",ndim_cell,dims_cell,ndim,dims); 
        par->m = mymex::set_field_cell(sol_struct,"m",ndim_cell,dims_cell,ndim,dims);

            delete[] dims_cell;
            delete[] ndim;
            delete[] dims[0];
            delete[] dims;

        logs::solve(1," output allocated\n");            

    } else if(type == 2){

        par->c = new double*[par->T];
        par->m = new double*[par->T];      

        for(size_t t = 0; t < par->T; t++){
            par->c[t] = (double*) mxGetPr(mxGetCell(mxGetField(prhs[1],0,"c"),t));
            par->m[t] = (double*) mxGetPr(mxGetCell(mxGetField(prhs[1],0,"m"),t));
        }

        logs::simulate(1," solution loaded\n");

    }

    ///////////////////////////
    // 4. outputs - simulate //
    ///////////////////////////

    if(type == 2){

        // a. struct
        const char *field_names[] = {"p", "phat",
                                     "xi", "psi", "xihat", "etahat",
                                     "A", 
                                     "MPC", "MPCP", "logC", "logY"};
        
        size_t num_fields = sizeof(field_names)/sizeof(*field_names);
        plhs[0] = mxCreateStructMatrix(1, 1, num_fields, field_names);
        auto sim_struct = plhs[0]; 

        // b. dimensions
        size_t ndim = 2;
        auto dims = new size_t[2];

            dims[0] = par->simN;
            dims[1] = par->simT;

        // c. elements
        par->sim_p      = mymex::set_field_double(sim_struct,"p",ndim,dims);
        par->sim_phat   = mymex::set_field_double(sim_struct,"phat",ndim,dims);
        par->sim_xi     = mymex::set_field_double(sim_struct,"xi",ndim,dims);
        par->sim_psi    = mymex::set_field_double(sim_struct,"psi",ndim,dims);        
        par->sim_xihat  = mymex::set_field_double(sim_struct,"xihat",ndim,dims);
        par->sim_etahat = mymex::set_field_double(sim_struct,"etahat",ndim,dims);
        par->sim_A      = mymex::set_field_double(sim_struct,"A",ndim,dims);
        par->sim_MPC    = mymex::set_field_double(sim_struct,"MPC",ndim,dims);
        par->sim_MPCP   = mymex::set_field_double(sim_struct,"MPCP",ndim,dims);
        par->sim_logC   = mymex::set_field_double(sim_struct,"logC",ndim,dims);
        par->sim_logY   = mymex::set_field_double(sim_struct,"logY",ndim,dims);

            delete[] dims;

    }


    ///////////////////////////
    // 4. outputs - logdiffs //
    ///////////////////////////

    if(type == 3){

        // a. struct
        const char *field_names[] = {"DlogY",
                                     "DlogY_lead1",
                                     "DlogY_lead2", 
                                     "DlogY_lead3", 
                                     "DlogY_lead4", 
                                     "DlogY_lead5",
                                     "DlogY_perm",
                                     "DlogY_lag1",
                                     "DlogY_lag2", 
                                     "DlogY_lag3", 
                                     "DlogY_lag4", 
                                     "DlogY_lag5",
                                     "DlogC",
                                     "DlogC_lead1",
                                     "DlogC_lead2",
                                     "DlogC_lead3",
                                     "DlogC_lead4",
                                     "DlogC_lead5"};
        
        size_t num_fields = sizeof(field_names)/sizeof(*field_names);
        plhs[0] = mxCreateStructMatrix(1, 1, num_fields, field_names);
        auto sim_struct = plhs[0]; 

        // b. dimensions
        size_t ndim = 2;
        auto dims = new size_t[2];

            dims[0] = par->simN;
            dims[1] = par->simT;

        // c. elements
        par->DlogY          = mymex::set_field_double(sim_struct,"DlogY",ndim,dims);
        par->DlogY_lead1    = mymex::set_field_double(sim_struct,"DlogY_lead1",ndim,dims);
        par->DlogY_lead2    = mymex::set_field_double(sim_struct,"DlogY_lead2",ndim,dims);
        par->DlogY_lead3    = mymex::set_field_double(sim_struct,"DlogY_lead3",ndim,dims);
        par->DlogY_lead4    = mymex::set_field_double(sim_struct,"DlogY_lead4",ndim,dims);
        par->DlogY_lead5    = mymex::set_field_double(sim_struct,"DlogY_lead5",ndim,dims);
        par->DlogY_perm     = mymex::set_field_double(sim_struct,"DlogY_perm",ndim,dims);
        par->DlogY_lag1     = mymex::set_field_double(sim_struct,"DlogY_lag1",ndim,dims);
        par->DlogY_lag2     = mymex::set_field_double(sim_struct,"DlogY_lag2",ndim,dims);
        par->DlogY_lag3     = mymex::set_field_double(sim_struct,"DlogY_lag3",ndim,dims);
        par->DlogY_lag4     = mymex::set_field_double(sim_struct,"DlogY_lag4",ndim,dims);
        par->DlogY_lag5     = mymex::set_field_double(sim_struct,"DlogY_lag5",ndim,dims);
        par->DlogC          = mymex::set_field_double(sim_struct,"DlogC",ndim,dims);
        par->DlogC_lead1    = mymex::set_field_double(sim_struct,"DlogC_lead1",ndim,dims);
        par->DlogC_lead2    = mymex::set_field_double(sim_struct,"DlogC_lead2",ndim,dims);
        par->DlogC_lead3    = mymex::set_field_double(sim_struct,"DlogC_lead3",ndim,dims);
        par->DlogC_lead4    = mymex::set_field_double(sim_struct,"DlogC_lead4",ndim,dims);
        par->DlogC_lead5    = mymex::set_field_double(sim_struct,"DlogC_lead5",ndim,dims);

            delete[] dims;        
        
    }

}


////////////////
// 3. destroy //
////////////////

void destroy(par_struct *par, size_t type){

    if(type == 1){
        delete[] par->K;
        if(par->Nxihat != 1){
            delete[] par->VsqrtD;
        }
    }
    if(type == 2){
        delete[] par->c;
        delete[] par->m;
    }

    #if ONEDIM == 1 
    delete[] par->Q_vec;
    #endif    

}

} // namespace