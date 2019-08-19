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

// c. solve
#include "includes\interpolate.cpp"     // setup and interpolant gatewys


////////////////
// 4. gateway //
////////////////

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

        logs::calc_logdiffs(0,"begining\n"); // clear log

    //////////////////
    // 5. setup par //
    //////////////////
    
    par_struct* par = new par_struct;
    par::setup(par,plhs,prhs,3);

    //////////////////////
    // 6. calc parallel //
    //////////////////////
    
    #pragma omp parallel num_threads(THREADS)
    {
	   
        // a. DlogY and DlogC    
        #pragma omp for      
        for(size_t t = 0; t < par->simT; t++){

            for(size_t i = 0; i < par->simN; i++){

                size_t now = t*par->simN + i;   
                size_t lag = (t-1)*par->simN + i;   

                par->DlogY[now] = t > 0 ? par->sim_logY[now]-par->sim_logY[lag] : NAN;
                par->DlogC[now] = t > 0 ? par->sim_logC[now]-par->sim_logC[lag] : NAN;                        

            }


        }

        // b. DlogY and DlogC lags and leads
        #pragma omp for
        for(size_t t = 0; t < par->simT; t++){

            for(size_t i = 0; i < par->simN; i++){
            
                size_t now   = t*par->simN + i;  

                size_t lag_5 = (t-5)*par->simN + i;
                size_t lag_4 = (t-4)*par->simN + i;
                size_t lag_3 = (t-3)*par->simN + i;
                size_t lag_2 = (t-2)*par->simN + i;
                size_t lag_1 = (t-1)*par->simN + i;  

                size_t lead_1 = (t+1)*par->simN + i;   
                size_t lead_2 = (t+2)*par->simN + i;
                size_t lead_3 = (t+3)*par->simN + i;
                size_t lead_4 = (t+4)*par->simN + i;
                size_t lead_5 = (t+5)*par->simN + i;

                par->DlogY_lag5[now] = t > 5 ? par->DlogY[lag_5] : NAN;
                par->DlogY_lag4[now] = t > 4 ? par->DlogY[lag_4] : NAN;
                par->DlogY_lag3[now] = t > 3 ? par->DlogY[lag_3] : NAN;
                par->DlogY_lag2[now] = t > 2 ? par->DlogY[lag_2] : NAN;
                par->DlogY_lag1[now] = t > 1 ? par->DlogY[lag_1] : NAN;

                par->DlogY_lead1[now] = t < par->simT-1-1 ? par->DlogY[lead_1] : NAN;
                par->DlogY_lead2[now] = t < par->simT-1-2 ? par->DlogY[lead_2] : NAN;
                par->DlogY_lead3[now] = t < par->simT-1-3 ? par->DlogY[lead_3] : NAN;
                par->DlogY_lead4[now] = t < par->simT-1-4 ? par->DlogY[lead_4] : NAN;
                par->DlogY_lead5[now] = t < par->simT-1-5 ? par->DlogY[lead_5] : NAN;

                par->DlogY_perm[now] = par->DlogY_lag1[now] + par->DlogY[now] + par->DlogY_lead1[now]; 

                par->DlogC_lead1[now] = t < par->simT-1-1 ? par->DlogC[lead_1] : NAN;
                par->DlogC_lead2[now] = t < par->simT-1-2 ? par->DlogC[lead_2] : NAN;
                par->DlogC_lead3[now] = t < par->simT-1-3 ? par->DlogC[lead_3] : NAN;
                par->DlogC_lead4[now] = t < par->simT-1-4 ? par->DlogC[lead_4] : NAN;
                par->DlogC_lead5[now] = t < par->simT-1-5 ? par->DlogC[lead_5] : NAN;                

            }

        }


    } // parallel


    /////////////////
    // 7. clean up //
    /////////////////

    par::destroy(par,3);
    delete par;

        logs::calc_logdiffs(1,"Done.\n");

        // clean assertions file
        FILE* log_file = fopen("log_assert.txt","w");
        fprintf(log_file,"\n");
        fclose(log_file);

} // simulate