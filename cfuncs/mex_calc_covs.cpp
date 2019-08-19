#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <string>
#include <omp.h>
#include "vectorclass/vectorclass.h"
#include "mex.h"
#include "matrix.h"

#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))
#define BOUND(X,A,B) MIN(MAX(X,A),B)
//#define THREADS MIN(MAXTHREADS, omp_get_max_threads()-1)
#define THREADS MAXTHREADS

#include "includes\HighResTimer_class.hpp" // timer class
#include "includes\assert.cpp"             // assert() function
#include "includes\logs.cpp"               // log:: functions
#include "includes\mymex.cpp"              // functions to interact with mex

double cov(double *A, double *B, double *g, size_t N, size_t do_full_output);


/////////////
// GATEWAY //
/////////////

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

	logs::calc_covs(0,"begining\n"); // clear log

	///////////////
	// 1. inputs //
	///////////////

	size_t Ntotal = mxGetNumberOfElements(mxGetField(prhs[0],0,"DlogY"));

	// a. income differences
    auto DlogY          = (double*) mxGetPr(mxGetField(prhs[0],0,"DlogY"));

	auto DlogY_lead1    = (double*) mxGetPr(mxGetField(prhs[0],0,"DlogY_lead1"));
	auto DlogY_lead2    = (double*) mxGetPr(mxGetField(prhs[0],0,"DlogY_lead2")); 
    auto DlogY_lead3    = (double*) mxGetPr(mxGetField(prhs[0],0,"DlogY_lead3")); 
    auto DlogY_lead4    = (double*) mxGetPr(mxGetField(prhs[0],0,"DlogY_lead4")); 
    auto DlogY_lead5    = (double*) mxGetPr(mxGetField(prhs[0],0,"DlogY_lead5"));

	auto DlogY_perm    = (double*) mxGetPr(mxGetField(prhs[0],0,"DlogY_perm"));

	auto DlogY_lag1    = (double*) mxGetPr(mxGetField(prhs[0],0,"DlogY_lag1"));
	auto DlogY_lag2    = (double*) mxGetPr(mxGetField(prhs[0],0,"DlogY_lag2")); 
    auto DlogY_lag3    = (double*) mxGetPr(mxGetField(prhs[0],0,"DlogY_lag3")); 
    auto DlogY_lag4    = (double*) mxGetPr(mxGetField(prhs[0],0,"DlogY_lag4")); 
    auto DlogY_lag5    = (double*) mxGetPr(mxGetField(prhs[0],0,"DlogY_lag5"));
    
    // b. consumption differences
    auto DlogC          = (double*) mxGetPr(mxGetField(prhs[0],0,"DlogC"));

	auto DlogC_lead1    = (double*) mxGetPr(mxGetField(prhs[0],0,"DlogC_lead1"));
	auto DlogC_lead2    = (double*) mxGetPr(mxGetField(prhs[0],0,"DlogC_lead2"));
	auto DlogC_lead3    = (double*) mxGetPr(mxGetField(prhs[0],0,"DlogC_lead3"));
	auto DlogC_lead4    = (double*) mxGetPr(mxGetField(prhs[0],0,"DlogC_lead4"));
	auto DlogC_lead5    = (double*) mxGetPr(mxGetField(prhs[0],0,"DlogC_lead5"));

	// c. full output
	size_t do_full_output = mxGetScalar(prhs[1]);

		// input assertions
		assert(0,Ntotal%4 == 0,"modulo(N,4) != 0 is not allowed");
	
		logs::calc_covs(1,"input loaded\n");


	///////////////
	// 2. output //
	///////////////

    // a. struct
	const char *field_names[] = {"yy",
								"yy_lead1",
								"yy_lead2",
								"yy_lead3",
								"yy_lead4",
								"yy_lead5",
								"cc",
								"cc_lead1",
								"cc_lead2",
								"cc_lead3",
								"cc_lead4",
								"cc_lead5",
								"cy_lag1",
								"cy_lag2",
								"cy_lag3",
								"cy_lag4",
								"cy_lag5",							 
								"cy",
								"cy_lead1",
								"cy_lead2",
								"cy_lead3",
								"cy_lead4",
								"cy_lead5",
								"yy_perm",
							    "cy_perm",
        							
        						"yy_vec",
								"yy_lead1_vec",
								"yy_lead2_vec",
								"yy_lead3_vec",
								"yy_lead4_vec",
								"yy_lead5_vec",
								"cc_vec",
								"cc_lead1_vec",
								"cc_lead2_vec",
								"cc_lead3_vec",
								"cc_lead4_vec",
								"cc_lead5_vec",
								"cy_lag1_vec",
								"cy_lag2_vec",
								"cy_lag3_vec",
								"cy_lag4_vec",
								"cy_lag5_vec",	
								"cy_vec",
								"cy_lead1_vec",
								"cy_lead2_vec",
								"cy_lead3_vec",
								"cy_lead4_vec",
								"cy_lead5_vec",
								"yy_perm_vec",
								"cy_perm_vec"};

		size_t num_fields = sizeof(field_names)/sizeof(*field_names);
	    plhs[0] = mxCreateStructMatrix(1, 1, num_fields, field_names);
	    auto output_struct = plhs[0];

	// b. dimensions	    
    size_t ndim = 1; 
    size_t dims[1];
    dims[0] = 1;

    // c. elements

	    // income
		auto yy		   = mymex::set_field_double(output_struct,"yy",ndim,dims);
		auto yy_lead1  = mymex::set_field_double(output_struct,"yy_lead1",ndim,dims);
		auto yy_lead2  = mymex::set_field_double(output_struct,"yy_lead2",ndim,dims);
		auto yy_lead3  = mymex::set_field_double(output_struct,"yy_lead3",ndim,dims);
		auto yy_lead4  = mymex::set_field_double(output_struct,"yy_lead4",ndim,dims);
		auto yy_lead5  = mymex::set_field_double(output_struct,"yy_lead5",ndim,dims);

		// consumption
		auto cc		   = mymex::set_field_double(output_struct,"cc",ndim,dims);
		auto cc_lead1  = mymex::set_field_double(output_struct,"cc_lead1",ndim,dims);
		auto cc_lead2  = mymex::set_field_double(output_struct,"cc_lead2",ndim,dims);
		auto cc_lead3  = mymex::set_field_double(output_struct,"cc_lead3",ndim,dims);
		auto cc_lead4  = mymex::set_field_double(output_struct,"cc_lead4",ndim,dims);
		auto cc_lead5  = mymex::set_field_double(output_struct,"cc_lead5",ndim,dims);

		// mixed
		auto cy_lag1   = mymex::set_field_double(output_struct,"cy_lag1",ndim,dims);
		auto cy_lag2   = mymex::set_field_double(output_struct,"cy_lag2",ndim,dims);
		auto cy_lag3   = mymex::set_field_double(output_struct,"cy_lag3",ndim,dims);
		auto cy_lag4   = mymex::set_field_double(output_struct,"cy_lag4",ndim,dims);
		auto cy_lag5   = mymex::set_field_double(output_struct,"cy_lag5",ndim,dims);
		auto cy		   = mymex::set_field_double(output_struct,"cy",ndim,dims);
		auto cy_lead1  = mymex::set_field_double(output_struct,"cy_lead1",ndim,dims);
		auto cy_lead2  = mymex::set_field_double(output_struct,"cy_lead2",ndim,dims);
		auto cy_lead3  = mymex::set_field_double(output_struct,"cy_lead3",ndim,dims);
		auto cy_lead4  = mymex::set_field_double(output_struct,"cy_lead4",ndim,dims);
		auto cy_lead5  = mymex::set_field_double(output_struct,"cy_lead5",ndim,dims);

		// permanent
		auto yy_perm   = mymex::set_field_double(output_struct,"yy_perm",ndim,dims);
		auto cy_perm   = mymex::set_field_double(output_struct,"cy_perm",ndim,dims);

	// d. fulle output elements
	if(do_full_output){
    	dims[0] = Ntotal;		
    }

	    // income
		auto yy_vec		   = mymex::set_field_double(output_struct,"yy_vec",ndim,dims);
		auto yy_lead1_vec  = mymex::set_field_double(output_struct,"yy_lead1_vec",ndim,dims);
		auto yy_lead2_vec  = mymex::set_field_double(output_struct,"yy_lead2_vec",ndim,dims);
		auto yy_lead3_vec  = mymex::set_field_double(output_struct,"yy_lead3_vec",ndim,dims);
		auto yy_lead4_vec  = mymex::set_field_double(output_struct,"yy_lead4_vec",ndim,dims);
		auto yy_lead5_vec  = mymex::set_field_double(output_struct,"yy_lead5_vec",ndim,dims);

		// consumption
		auto cc_vec		   = mymex::set_field_double(output_struct,"cc_vec",ndim,dims);
		auto cc_lead1_vec  = mymex::set_field_double(output_struct,"cc_lead1_vec",ndim,dims);
		auto cc_lead2_vec  = mymex::set_field_double(output_struct,"cc_lead2_vec",ndim,dims);
		auto cc_lead3_vec  = mymex::set_field_double(output_struct,"cc_lead3_vec",ndim,dims);
		auto cc_lead4_vec  = mymex::set_field_double(output_struct,"cc_lead4_vec",ndim,dims);
		auto cc_lead5_vec  = mymex::set_field_double(output_struct,"cc_lead5_vec",ndim,dims);

		// mixed
		auto cy_lag1_vec   = mymex::set_field_double(output_struct,"cy_lag1_vec",ndim,dims);
		auto cy_lag2_vec   = mymex::set_field_double(output_struct,"cy_lag2_vec",ndim,dims);
		auto cy_lag3_vec   = mymex::set_field_double(output_struct,"cy_lag3_vec",ndim,dims);
		auto cy_lag4_vec   = mymex::set_field_double(output_struct,"cy_lag4_vec",ndim,dims);
		auto cy_lag5_vec   = mymex::set_field_double(output_struct,"cy_lag5_vec",ndim,dims);
		auto cy_vec   	   = mymex::set_field_double(output_struct,"cy_vec",ndim,dims);
		auto cy_lead1_vec  = mymex::set_field_double(output_struct,"cy_lead1_vec",ndim,dims);
		auto cy_lead2_vec  = mymex::set_field_double(output_struct,"cy_lead2_vec",ndim,dims);
		auto cy_lead3_vec  = mymex::set_field_double(output_struct,"cy_lead3_vec",ndim,dims);
		auto cy_lead4_vec  = mymex::set_field_double(output_struct,"cy_lead4_vec",ndim,dims);
		auto cy_lead5_vec  = mymex::set_field_double(output_struct,"cy_lead5_vec",ndim,dims);

		// permanent
		auto yy_perm_vec   = mymex::set_field_double(output_struct,"yy_perm_vec",ndim,dims);
		auto cy_perm_vec   = mymex::set_field_double(output_struct,"cy_perm_vec",ndim,dims);

		logs::calc_covs(1,"output allocated\n");

	/////////////
	// 3. covs //
	/////////////

    #pragma omp parallel num_threads(THREADS)
    {

		#pragma omp single
		{
			
			////////////
			// income //
			////////////

	    	#pragma omp task
	    	yy[0] = cov(DlogY, DlogY, yy_vec, Ntotal, do_full_output);


			#pragma omp task
	    	yy_lead1[0] = cov(DlogY, DlogY_lead1, yy_lead1_vec, Ntotal, do_full_output); 
			
			#pragma omp task
			yy_lead2[0] = cov(DlogY, DlogY_lead2, yy_lead2_vec, Ntotal, do_full_output);

			#pragma omp task
			yy_lead3[0] = cov(DlogY, DlogY_lead3, yy_lead3_vec, Ntotal, do_full_output);
			
			#pragma omp task
			yy_lead4[0] = cov(DlogY, DlogY_lead4, yy_lead4_vec, Ntotal, do_full_output);
			
			#pragma omp task
			yy_lead5[0] = cov(DlogY, DlogY_lead5, yy_lead5_vec, Ntotal, do_full_output);
			

			/////////////////
			// consumption //
			/////////////////
			
			#pragma omp task
			cc[0] = cov(DlogC, DlogC, cc_vec, Ntotal, do_full_output);

					 
			#pragma omp task
			cc_lead1[0] = cov(DlogC, DlogC_lead1, cc_lead1_vec, Ntotal, do_full_output);
			
			#pragma omp task
			cc_lead2[0] = cov(DlogC, DlogC_lead2, cc_lead2_vec, Ntotal, do_full_output);
			
			#pragma omp task
			cc_lead3[0] = cov(DlogC, DlogC_lead3, cc_lead3_vec, Ntotal, do_full_output);
			
			#pragma omp task
			cc_lead4[0] = cov(DlogC, DlogC_lead4, cc_lead4_vec, Ntotal, do_full_output);
			
			#pragma omp task
			cc_lead5[0] = cov(DlogC, DlogC_lead5, cc_lead5_vec, Ntotal, do_full_output);
			

			#pragma omp task
			cy_lag1[0] = cov(DlogC, DlogY_lag1, cy_lag1_vec, Ntotal, do_full_output);
	
			#pragma omp task
			cy_lag2[0] = cov(DlogC, DlogY_lag2, cy_lag2_vec, Ntotal, do_full_output);

			#pragma omp task
			cy_lag3[0] = cov(DlogC, DlogY_lag3, cy_lag3_vec, Ntotal, do_full_output);

			#pragma omp task
			cy_lag4[0] = cov(DlogC, DlogY_lag4, cy_lag4_vec, Ntotal, do_full_output);
	
			#pragma omp task
			cy_lag5[0] = cov(DlogC, DlogY_lag5, cy_lag5_vec, Ntotal, do_full_output);
			 

			///////////
			// mixed //
			///////////
			
			#pragma omp task
			cy[0] = cov(DlogC, DlogY, cy_vec, Ntotal, do_full_output);
				

			#pragma omp task
			cy_lead1[0] = cov(DlogC, DlogY_lead1, cy_lead1_vec, Ntotal, do_full_output);
			
			#pragma omp task
			cy_lead2[0] = cov(DlogC, DlogY_lead2, cy_lead2_vec, Ntotal, do_full_output);
			
			#pragma omp task
			cy_lead3[0] = cov(DlogC, DlogY_lead3, cy_lead3_vec, Ntotal, do_full_output);
			
			#pragma omp task
			cy_lead4[0] = cov(DlogC, DlogY_lead4, cy_lead4_vec, Ntotal, do_full_output);
			
			#pragma omp task
			cy_lead5[0] = cov(DlogC, DlogY_lead5, cy_lead5_vec, Ntotal, do_full_output);
			         	

			///////////////
			// permanent //
			///////////////

			#pragma omp task
	    	yy_perm[0] = cov(DlogY, DlogY_perm, yy_perm_vec, Ntotal, do_full_output);

			#pragma omp task
	    	cy_perm[0] = cov(DlogC, DlogY_perm, cy_perm_vec, Ntotal, do_full_output);


		}

	}

		logs::calc_covs(1,"Done.\n");

}

//////////
// MAIN //
//////////

double cov(double *A, double *B, double *g, size_t N, size_t do_full_output){

	// a. determine sums and sample
	double mean_A = 0.0;
	double mean_B = 0.0;
	size_t Nactive = 0;
	for(size_t i = 0; i < N; i+=4){
		
		// a. load 
		Vec4d A_vec, B_vec;
		A_vec.load(&A[i]);
		B_vec.load(&B[i]);

		// b. not nan
		Vec4db I_A_vec = !is_nan(A_vec);
		Vec4db I_B_vec = !is_nan(B_vec);
		Vec4db I_vec = I_A_vec & I_B_vec;

		// c. replace nans with 0
		Vec4d A_selected = select(I_vec, A_vec, 0);
		Vec4d B_selected = select(I_vec, B_vec, 0);

		// c. sum
		mean_A += horizontal_add(A_selected);
		mean_B += horizontal_add(B_selected);
		Nactive += horizontal_count(I_vec);

	}

	// b. means
	mean_A /= (double)Nactive;
	mean_B /= (double)Nactive;

	// c. covariance
	double cov = 0.0;
	for(size_t i = 0; i < N; i += 4){
		
		// a. load 
		Vec4d A_vec, B_vec;
		A_vec.load(&A[i]);
		B_vec.load(&B[i]);

		// b. not nan
		Vec4db I_A_vec = !is_nan(A_vec);
		Vec4db I_B_vec = !is_nan(B_vec);
		Vec4db I_vec = I_A_vec & I_B_vec;

		// b. calculate
		Vec4d AB_vec = (A_vec-mean_A)*(B_vec-mean_B);

		// c. select
		Vec4d AB_selected = select(I_vec, AB_vec, 0);

		// d. sum
		cov += horizontal_add(AB_selected);

		// e. save
		if(do_full_output){
			Vec4d g_vec = select(I_vec, AB_vec, NAN);
			g_vec.store(&g[i]);
		}
		
	}

	// d. result
	if(Nactive == 0){
		return NAN;
	} else {
		return cov/(double)(Nactive);
 	}
 	
}