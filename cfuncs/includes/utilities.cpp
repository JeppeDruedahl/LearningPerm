namespace utilities {

////////////////
// 1. utility //
////////////////

double marg_u(double C, par_struct *par)
{
    return pow(C,-par->rho);
    
}

double inv_marg_u(double u, par_struct *par)
{
    return pow(u,-1.0/par->rho);
}

}