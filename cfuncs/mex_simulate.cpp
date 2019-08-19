//////////////////////////
// 1. external includes //
//////////////////////////

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <string>
#include <omp.h>
#include "mex.h"
#include "matrix.h"


//////////////////////////
// 2. define statements //
//////////////////////////

#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))
#define BOUND(X,A,B) MIN(MAX(X,A),B)
#define THREADS MAXTHREADS
#define ONEDIM 0

//////////////////////////
// 3. internal includes //
//////////////////////////

// a. generic
#include "includes\HighResTimer_class.hpp" // timer class
#include "includes\assert.cpp"             // assert() function
#include "includes\logs.cpp"               // log:: functions
#include "includes\linear_interp.cpp"      // linear_interp:: functions
#include "includes\index.cpp"              // index:: functions
#include "includes\mymex.cpp"              // functions to interact with mex

// b. basic
#include "includes\par_struct.cpp"  // define par_struct + setup/destroy functions
#include "includes\sol_struct.cpp"  // define sol_struct + setup/destroy functions

// c. solve
#include "includes\interpolate.cpp"     // setup and interpolant gatewys


////////////////
// 4. gateway //
////////////////

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

    HighResTimer timer_all, timer;
    timer_all.StartTimer();

        logs::simulate(0,"Starting simulation...\n");


    //////////////////
    // 5. setup par //
    //////////////////
    
    par_struct* par = new par_struct;
    par::setup(par,plhs,prhs,2);

        logs::simulate(1,"Setup completed.\n");

        // easy access to outputs
        auto p       = par->sim_p;
        auto phat    = par->sim_phat;    
        auto xi      = par->sim_xi;           
        auto psi     = par->sim_psi;                   
        auto xihat   = par->sim_xihat;           
        auto etahat  = par->sim_etahat;
        auto A       = par->sim_A;
        auto MPC     = par->sim_MPC;
        auto MPCP    = par->sim_MPCP;
        auto logC    = par->sim_logC;
        auto logY    = par->sim_logY;

    /////////////////
    // 6. parallel //
    /////////////////
        
    #pragma omp parallel num_threads(THREADS)
    {

        auto sol = new sol_struct;
        sol::setup(par,sol);  


    //////////////////
    // 7. time loop //
    //////////////////

    for(size_t t = 0; t < par->simT; t++){ // forward through time

        #pragma omp master
        logs::simulate(2,"t = %d\n",t);
        
        #pragma omp master
        timer.StartTimer();

        // create interpolants
        if(t == 0){
            interpolate::setup(t,par,sol,1);
        } else {
            interpolate::update(t,par,sol);
        }

    /////////////////////////
    // 8. individuals loop //
    /////////////////////////

        logs::simulate(3,"individuals loop\n");

    #pragma omp for schedule(guided)
    for(size_t i = 0; i < par->simN; i++){
	
            logs::simulate(4,"i = %d\n",i);

        // a. index
        size_t index = t*par->simN + i;    
	    size_t index_lag = (t-1)*par->simN + i;

        // b. lag
        double p_lag, phat_lag, A_lag, xi_lag, xihat_lag;
        if(t == 0){
            
            p_lag    = 0.0;
            phat_lag = p_lag + sqrt(par->q0)*par->draws_q0[i];
            A_lag    = par->A0;//0.0;

            xi_lag    = 0.0;
            xihat_lag = 0.0;

        } else {
            
            p_lag    = p[index_lag];
            phat_lag = phat[index_lag];
            A_lag    = A[index_lag];
        
            xi_lag    = xi[index_lag];
            xihat_lag = xihat[index_lag];

        }    

        // c. income
        double y,z;
        if(t < par->TR){

            psi[index] = par->mu_psi + par->sigma_psi*par->draws_psi[index];
            xi[index] = par->mu_xi + par->sigma_xi*par->draws_xi[index];
            double eps = 0.0;
            if(mxIsInf(par->sigma_eps) == 0){
                eps = par->sigma_eps*par->draws_eps[index];
            }
            p[index] = par->logG[t] + par->alpha*p_lag + psi[index];
            y = p[index] + xi[index] + par->omega*xi_lag;
            z = p[index] + eps;

        } else {
            
            p[index] = p_lag;
            y = p[index] + par->logkappa;
            z = p[index];
        
        }
        double Y = exp(y); 

        // c. belief 
        if(t < par->TR){
            
            double y_predict;
            if(par->Nxihat == 1){
            
                double p_predict = par->logG[t] + par->alpha*phat_lag + par->mu_psi;
                y_predict = p_predict + par->mu_xi;
                double z_predict = p_predict;

                phat[index] = p_predict + par->K[t][0]*(y-y_predict) + par->K[t][1]*(z-z_predict);
            
            } else {

                double p_predict   = par->logG[t] + par->alpha*phat_lag + par->mu_psi;
                double eta_predict = par->omega*xihat_lag + par->mu_xi;
                double xi_predict  = par->mu_xi;
                y_predict   = p_predict + eta_predict;
                double z_predict   = p_predict;
                
                double Delta_1 = y-y_predict;
                double Delta_2 = z-z_predict;

                phat[index]  = p_predict + par->K[t][0]*Delta_1 + par->K[t][3]*Delta_2;
                etahat[index]  = eta_predict + par->K[t][1]*Delta_1 + par->K[t][4]*Delta_2;                
                xihat[index] = xi_predict + par->K[t][2]*Delta_1 + par->K[t][5]*Delta_2;

            }

        } else {

            phat[index] = p[index];            
        
        }
        double Phat = exp(phat[index]);

        // d. states
        double M = par->R*A_lag + Y;
        double m = M/Phat; 

        // e. consumption
        double c = interpolate::evaluate_single(par,sol,m,phat[index],xihat[index]);             
        double C = c*Phat;
        
        // f. end-of-period assets
        A[index] = M-C;        
        
        // g. MPC and MPCP 
        if(par->do_MPC_MPCP){   
            
            // i. MPC
            double M_add = M + 0.01;

            double m_add = M_add/Phat;
            double c_add = interpolate::evaluate_single(par,sol,m_add,phat[index],xihat[index]); 
            double C_add = c_add*Phat;

            MPC[index] = (C_add - C) / (M_add - M);

            // ii. MPCP
            M_add = M + 0.01;
            double Phat_add  = Phat + 0.01;

            m_add = M_add/Phat_add;
            c_add = interpolate::evaluate_single(par,sol,m_add,phat[index],xihat[index]); 
            C_add = c_add*Phat_add;

            MPCP[index] = (C_add - C) / (M_add - M);
            
        }
        
        // 4. measurement error
        logC[index] = log(C) + par->sigma_eta_c*par->draws_eta_c[index];
        logY[index] = y + par->sigma_eta_y*par->draws_eta_y[index];


    } // N
    } // T        

        interpolate::destroy(par,sol);
        
        sol::destroy(sol);
        delete sol;
	
    } // parallel


    //////////////////
    // 9. clean up //
    //////////////////

    par::destroy(par,2);
    delete par;

        double time_all = timer_all.StopTimer();
        logs::simulate(1,"Time: %5.2f secs\n",time_all);
        logs::simulate(1,"Done.\n");

        // clean assertions file
        FILE* log_file = fopen("log_assert.txt","w");
        fprintf(log_file,"\n");
        fclose(log_file);

} // simulate