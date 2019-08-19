namespace linear_interp {

typedef struct set_struct // set is short for settings
{
    // a. input
    size_t dimx;    // number of dimensions in x
    size_t dimy;    // number of dimensions in y
    size_t *Nx;     // lenght of each x-dimension

    double **x;  // pointers to grid vectors
    double **y;  // pointers to value arrays

    size_t Nxi;     // number of points to be evaluated

    size_t ncube;   // number of corners in hypercube of dim = dimx-1

    // b. containers - standard
    size_t *facs, *add, *pos_left, *indexes;
    double *reldiff, *weights;
    double **cs;

    // c. containers - vectorized
    size_t *pos_left_vec;
    double *reldiff_vec;

} set_struct;

void create(linear_interp::set_struct *set, 
            size_t dimx, size_t dimy, size_t *Nx, 
            double **x, double **y, 
            size_t Nxi)
{
    // a. input
    set->dimx  = dimx;
    set->dimy  = dimy;    
    set->Nxi   = Nxi;
    set->ncube = (size_t)pow(2.0,(double)(dimx-1));

    set->Nx = new size_t[dimx];
    set->x  = new double*[dimx];
    for(size_t j = 0; j < dimx; j++){
        set->Nx[j] = Nx[j];
        set->x[j]  = x[j];
    }

    set->y = new double*[dimy];
    for(size_t j = 0; j < dimy; j++){
        set->y[j] = y[j];
    }

    // b. containers - standard
    set->facs     = new size_t[dimx];
    set->add      = new size_t[set->ncube*dimx];
    set->pos_left = new size_t[dimx];
    set->indexes  = new size_t[set->ncube];
    
    set->reldiff  = new double[dimx];
    set->weights  = new double[set->ncube*2];

    // c. containers - vectorized
    set->pos_left_vec = new size_t[Nxi];
    set->reldiff_vec  = new double[Nxi];

    // d. position factors    
    for(size_t i = 0; i < dimx; i++){
        set->facs[i] = 1;
        for(size_t j = 0; j < i; j++){
            set->facs[i] *= Nx[j];
        }
    }

    // e. add (ordering of hypercube corners is important)
    for(size_t i = 0; i < set->ncube; i++){

        size_t *add = &set->add[i*dimx];

        if(i == 0){
            for(size_t j = 0; j < dimx; j++){ add[j] = 0;}
        } else {
            for(size_t j = 0; j < dimx; j++){ add[j] = set->add[(i-1)*dimx + j]; }
            for(size_t j = 1; j < dimx; j++){
                if(add[j] == 0){
                    add[j] = 1;
                    for(size_t k = 1; k < j; k++){ add[k] = 0; }
                    break;
                }
            }
        } // if i == 0
    } // i

    // f. dimension reduction
    set->cs = new double*[dimx];
    for(size_t j = 0; j < dimx; j++){
        set->cs[j] = new double[(size_t)pow(2.0,(double)(dimx-1-j))];
    }
}

void destroy(linear_interp::set_struct *set)
{

    // a. input
    delete[] set->Nx;    
    delete[] set->x;
    delete[] set->y; 

    // b. containers - standard
    delete[] set->facs;
    delete[] set->add;
    delete[] set->pos_left;
    delete[] set->indexes;

    delete[] set->reldiff;
    delete[] set->weights;

    // d. container - vectorized
    delete[] set->pos_left_vec;
    delete[] set->reldiff_vec;

    // f. dimension reduction
    for(size_t j = 0; j < set->dimx; j++){
        delete[] set->cs[j];
    }
    delete[] set->cs;

    // d. itself   
    delete set;
}

size_t binary_search(size_t imin, size_t Nx, double *x, double xi)
{
    size_t imid, half;

    // a. checks
    if(xi <= x[0]){
        return 0;
    } else if(xi >= x[Nx-2]) {
        return Nx-2;
    }

    // b. binary search
    while((half = Nx/2)){
        imid = imin + half;
        imin = (x[imid] <= xi) ? imid:imin;
        Nx  -= half;
    }

    return imin;
}

size_t find_index_left(linear_interp::set_struct *set, size_t *add)
{

    // a. unpack
    auto dimx     = set->dimx;
    auto facs     = set->facs;
    auto pos_left = set->pos_left;

    // b. calculate
    size_t index_left = 0;
    for(size_t j = 1; j < dimx; j++){
        index_left += facs[j]*(pos_left[j] + add[j]);
    }

    return index_left;
}

void find_reldiff(linear_interp::set_struct *set, double *xi)
{

    // a. unpack
    auto dimx     = set->dimx;
    auto x        = set->x;
    auto Nx       = set->Nx;
    auto pos_left = set->pos_left;
    auto reldiff  = set->reldiff;

    // b. find relative differences
    for(size_t d = 0; d < dimx; d++){

        // a. x-vector in dimension d
        double *xvec = x[d];

        // b. position to the left of xi in dimension d
        pos_left[d] = linear_interp::binary_search(0, Nx[d], xvec, xi[d]);

        // c. relative position if xi between neighboring points
        double denom = (xvec[pos_left[d]+1] - xvec[pos_left[d]]);
        reldiff[d] = (xi[d] - xvec[pos_left[d]]) / denom;

    }
}

void calc_val(linear_interp::set_struct *set, double *yi)
{

    // a. unpack
    auto dimy     = set->dimy;
    auto ncube    = set->ncube;
    auto add      = set->add;
    auto dimx     = set->dimx;    
    auto facs     = set->facs;
    auto pos_left = set->pos_left;
    auto indexes  = set->indexes;
    auto weights  = set->weights;
    auto reldiff  = set->reldiff;
    auto y        = set->y;
    
    // b. faster version
    if(set->Nxi == 1 && set->dimy == 1){

        // i. preliminaries
        size_t n = set->ncube;
        double *c = set->cs[0];

        // ii. function values at corner        
        size_t d = 0;
        for(size_t i = 0; i < n; i++){    

            size_t index_left = 0;
            for(size_t j = 0; j < dimx; j++){
                index_left += facs[j]*(pos_left[j] + add[i*dimx + j]);
            }

            c[i] = y[0][index_left]*(1.0-reldiff[d])+y[0][index_left+1]*reldiff[d];

        }     

        // iii. dimension reduction
        for(size_t d = 1; d < dimx; d++){

            // a. number of corners
            n /= 2;

            // b. calculate c in current dimension
            double *cprev = c;
            c     = set->cs[d];
            for(size_t i = 0; i < n; i++){
                c[i] = cprev[2*i]*(1.0-reldiff[d])+cprev[2*i+1]*reldiff[d];
            }

        }

        yi[0] = c[0];

    // c. might be a bit slower, but work for vectorized interpolants
    } else {

        // i. initialize at zero
        for(size_t l = 0; l < dimy; l++){
            yi[l] = 0;
        }

        // ii. calculate
        for(size_t i = 0; i < ncube; i++){

            indexes[i] = linear_interp::find_index_left(set, &add[i*dimx]);
            size_t index0 = facs[0]*(pos_left[0] + add[i*dimx+0]);

            weights[i*2 + 0] = 1.0;
            weights[i*2 + 1] = 1.0;
            for(size_t j = 1; j < dimx; j++){
                if(add[i*dimx + j] == 1){
                    weights[i*2 + 0] *= reldiff[j];
                    weights[i*2 + 1] *= reldiff[j];
                } else {
                    weights[i*2 + 0] *= (1.0-reldiff[j]);
                    weights[i*2 + 1] *= (1.0-reldiff[j]);
                }
            }

            for(size_t l = 0; l < dimy; l++){
                yi[l] += y[l][indexes[i]+index0]*(weights[i*2+0]*(1.0-reldiff[0]));
                yi[l] += y[l][indexes[i]+index0+1]*(weights[i*2+1]*reldiff[0]);
            }
        }

    }
}

void evaluate(linear_interp::set_struct *set,
              double *yi, double *xi, double *xi_vec)
{
    //////////////
    // 1. k = 0 //
    //////////////

    linear_interp::find_reldiff(set, xi);
    linear_interp::calc_val(set, yi);
    
    if(set->Nxi == 1){return;};

        // move result to vector format
        for(size_t l = 0; l < set->dimy; l++){
            yi[l*set->Nxi + 0] = yi[l];
        }

    //////////////
    // 2. k > 0 //
    //////////////

    // a. unpack
    auto x            = set->x;
    auto pos_left     = set->pos_left;
    auto pos_left_vec = set->pos_left_vec;
    auto Nxi          = set->Nxi;    
    auto dimy         = set->dimy;
    auto Nx           = set->Nx;    
    auto reldiff_vec  = set->reldiff_vec;
    auto indexes      = set->indexes;
    auto weights      = set->weights;
    auto ncube        = set->ncube;
    auto add          = set->add;
    auto dimx         = set->dimx;    
    auto y            = set->y;

    // b. loop through vector
    double *xvec = x[0];
    pos_left_vec[0] = pos_left[0];
    for(size_t k = 1; k < Nxi; k++){

        // i. initialize at zero
        for(size_t l = 0; l < dimy; l++){
            yi[l*Nxi + k] = 0.0;
        }

        // ii. position to the left of xi in dimension 0
        size_t i = pos_left_vec[k-1];
        while(xi_vec[k] > xvec[i+1] && i < Nx[0]-2){ i++;}
        pos_left_vec[k] = i;

        // ii. relative position if xi between neighboring points
        double denom = (xvec[pos_left_vec[k]+1] - xvec[pos_left_vec[k]]);
        reldiff_vec[k] = (xi_vec[k] - xvec[pos_left_vec[k]]) / denom;

    }

    // c. evaluate
    for(size_t l = 0; l < dimy; l++){
    for(size_t i = 0; i < ncube; i++){
        for(size_t k = 1; k < Nxi; k++){
            size_t index0 = (pos_left_vec[k] + add[i*dimx+0]);        
            yi[l*Nxi + k] += y[l][indexes[i]+index0]*(weights[i*2+0]*(1.0-reldiff_vec[k]));
            yi[l*Nxi + k] += y[l][indexes[i]+index0+1]*(weights[i*2+1]*reldiff_vec[k]);
        }
    }
    }

}

} // namespace