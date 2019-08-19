namespace interpolate {

//////////////
// 1. setup //
//////////////

void setup(size_t t, par_struct *par, sol_struct *sol, size_t Nxi){

	size_t i_cell = t;
	size_t i_grid = 0;

	// a. grids				
	size_t dimx;
	size_t Nx[3];
	double* x[3];

	if(par->Nphat == 1){
		dimx = 1;
		Nx[0] = par->Nm;
		x[0]  = &par->m[i_cell][i_grid];
	} else if(par->Nxihat == 1){
		dimx = 2;
		Nx[0] = par->Nm;
		Nx[1] = par->Nphat;
		x[0]  = &par->m[i_cell][i_grid];
		x[1]  = par->grid_phat;
	} else {
		dimx = 3;
		Nx[0] = par->Nm;
		Nx[1] = par->Nphat;
		Nx[2] = par->Nxihat;
		x[0]  = &par->m[i_cell][i_grid];
		x[1]  = par->grid_phat;
		x[2]  = par->grid_xihat;
	}
		
	// b. values
	size_t dimy;
	double* y[1];

	dimy = 1;
	y[0] = &par->c[i_cell][i_grid];

	// c. create		
	sol->interp = new linear_interp::set_struct;
	linear_interp::create(sol->interp, dimx, dimy, Nx, x, y, Nxi);

}

void setup_to_common(par_struct *par, sol_struct *sol){

	// a. grids				
	size_t dimx;
	size_t Nx[1];
	double* x[1];

	dimx = 1;
	Nx[0] = par->Na;
	x[0]  = sol->m_vec;
		
	// b. values
	size_t dimy;
	double* y[1];

	dimy = 1;
	y[0] = sol->c_vec;

	// c. create		
	sol->interp_to_common = new linear_interp::set_struct;
	linear_interp::create(sol->interp_to_common, dimx, dimy, Nx, x, y, par->Nm);

}

void update(size_t t, par_struct *par, sol_struct *sol){

	size_t i_cell = t;
	size_t i_grid = 0;

	sol->interp->x[0] = &par->m[i_cell][i_grid];   
	sol->interp->y[0] = &par->c[i_cell][i_grid];

}


void destroy(par_struct *par, sol_struct *sol){

	linear_interp::destroy(sol->interp);
	
}

void destroy_to_common(par_struct *par, sol_struct *sol){

	linear_interp::destroy(sol->interp_to_common);
	
}

/////////////////////
// 2. interpolants //
/////////////////////

void evaluate(par_struct *par, sol_struct *sol){

	double xi[3];
	xi[0] = sol->m_plus_vec[0];
	if(par->Nphat != 1){
		xi[1] = sol->phat_plus;
	} 
	if(par->Nxihat != 1){
		xi[2] = sol->xihat_plus;
	}
	linear_interp::evaluate(sol->interp,sol->c_plus_vec,xi,sol->m_plus_vec);

}

double evaluate_single(par_struct *par, sol_struct *sol, double m, double phat, double xihat){

	double xi[3];
	double c[1];

	xi[0] = m;
	if(par->Nphat != 1){
		xi[1] = phat;
	} 
	if(par->Nxihat != 1){
		xi[2] = xihat;
	}	
	linear_interp::evaluate(sol->interp,c,xi,nullptr);

	return c[0];

}

void evaluate_to_common(par_struct *par, sol_struct *sol, double* c){

	double xi[1];	
	xi[0] = par->grid_m[sol->t*par->Nm+0];
	linear_interp::evaluate(sol->interp_to_common,c,xi,&par->grid_m[sol->t*par->Nm+0]);

}

} // namespace