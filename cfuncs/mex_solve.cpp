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
#include "includes\utilities.cpp"   // transformation and utility functions

// c. solve
#include "includes\interpolate.cpp"     // setup and interpolant gatewys
#include "includes\EGM.cpp"             // EGM


////////////////
// 4. gateway //
////////////////

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
     
    HighResTimer timer_all, timer;
    timer_all.StartTimer();

        logs::solve(0,"Solving model.\n"); // reset log file
        

    //////////////////
    // 5. setup par //
    //////////////////
    
    par_struct* par = new par_struct;
    par::setup(par,plhs,prhs,1);

        logs::solve(1,"Setup completed.\n");

    /////////////////
    // 6. parallel //
    /////////////////

    #pragma omp parallel num_threads(THREADS)
    {

        // a. setup workers
        auto sol = new sol_struct;
        sol::setup(par,sol);

        // b. time loop
        for(size_t t = par->T; t-- > par->tmin;){

                #pragma omp master
                logs::solve(2,"t = %d ", t);

                #pragma omp master
                timer.StartTimer();

            sol->t = t;

            // i. interpolants
            if(t == par->T-2){
                interpolate::setup(t+1,par,sol,par->Na);
                if(par->Nphat != 1){
                   interpolate::setup_to_common(par,sol); 
                }
            } else if(t < par->T-2){
                interpolate::update(t+1,par,sol);                
            }

            // ii.EGM
            EGM::all(t,par,sol);
  
                #pragma omp master 
                {
                double time = timer.StopTimer();
                logs::solve(2,"time = %3.2f secs\n", time);                                                                
                }

        }

        // c. clean up workers
        interpolate::destroy(par,sol);        
        if(par->Nphat != 1){
            interpolate::destroy_to_common(par,sol); 
        }
        sol::destroy(sol);
        delete sol;

    } // parallel


    /////////////////
    // 7. clean up //
    /////////////////
    
    par::destroy(par,1);
    delete par;

        double time_all = timer_all.StopTimer();
        logs::solve(1,"Time: %5.2f secs\n",time_all);
        logs::solve(1,"Done.\n");

        // clean assertions file
        FILE* log_file = fopen("log_assert.txt","w");
        fprintf(log_file,"\n");
        fclose(log_file);

} // mex gateway