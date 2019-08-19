//////////////////
// 1. variables //
//////////////////

typedef struct
{
 
    // a. par
    par_struct *par;

    // b. interpolants
    linear_interp::set_struct *interp, *interp_to_common;

    // c. info
    size_t t;
    double phat, xihat, Phat;

    // d. EGM
    size_t i_shock;
    double *a_vec, *c_vec, *m_vec;
    double Psi_plus, phat_plus, xihat_plus, Phat_plus;
    double *Q_vec;
    double *m_plus_vec, *c_plus_vec;

} sol_struct;


namespace sol {

//////////////
// 2. setup //
//////////////

void setup(par_struct *par, sol_struct *sol){

    // a. par
    sol->par = par;      

    // c. EGM
    if(par->Nphat != 1){
        sol->m_vec = new double[par->Na];
        sol->c_vec = new double[par->Na];
    }
    sol->m_plus_vec = new double[par->Na];
    sol->c_plus_vec = new double[par->Na];

    #if ONEDIM == 0 
    sol->Q_vec = new double[par->Na];
    #else
    sol->Q_vec = par->Q_vec;    
    #endif        

}


////////////////
// 3. destroy //
////////////////

void destroy(sol_struct *sol){

    // c. EGM
    if(sol->par->Nphat != 1){
        delete[] sol->m_vec;
        delete[] sol->c_vec;        
    }
    delete[] sol->m_plus_vec;
    delete[] sol->c_plus_vec;

    #if ONEDIM == 0 
    delete[] sol->Q_vec;
    #endif        

}

} // namespace