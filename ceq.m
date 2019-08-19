classdef ceq
methods(Static)

function par = estimate(par,data)
    
    % a. point estimate
    
        yy       = data.moms.yy;
        cc_lead1 = data.moms.cc_lead1;
        yy_lead1 = data.moms.yy_lead1;
        yy_perm  = data.moms.yy_perm;
        cy_lead1 = data.moms.cy_lead1;

    par = ceq.point(par,yy,cc_lead1,yy_lead1,yy_perm,cy_lead1);
    
    % b. bootstrap
        
        % allocate
        for j = 1:numel(par.est_par)
            varname = sprintf('%s_bs',par.est_par{j});
            par.(varname) = nan(par.Nbootstraps_ceq,1);
        end    
    
    for i = 1:par.Nbootstraps_ceq
        
        par_bs = par;
        
        % i. moments
        yy       = data.moms.yy_bs(i);
        cc_lead1 = data.moms.cc_lead1_bs(i);
        yy_lead1 = data.moms.yy_lead1_bs(i);
        yy_perm  = data.moms.yy_perm_bs(i);
        cy_lead1 = data.moms.cy_lead1_bs(i);
        
        % ii. estimate
        par_bs = ceq.point(par_bs,yy,cc_lead1,yy_lead1,yy_perm,cy_lead1);
    
        % iii. save
        for j = 1:numel(par.est_par)
            varname = sprintf('%s_bs',par.est_par{j});
            par.(varname)(i) = par_bs.(par.est_par{j});
        end    
        
    end
        
    % c. standard errors    
    for j = 1:numel(par.est_par)
        varname_lhs = sprintf('%s_se',par.est_par{j});
        varname_rhs = sprintf('%s_bs',par.est_par{j});        
        par.(varname_lhs) = nanstd(par.(varname_rhs));
    end    
        
end
function par = point(par,yy,cc_lead1,yy_lead1,yy_perm,cy_lead1)
    
    % a. measurement error
    if ~isnan(par.meas_y_frac)
        par.sigma_eta_y = sqrt(0.5*par.meas_y_frac*yy); 
    end
    par.var_eta_c = max(-cc_lead1,0);

        % standard deviations
        par.sigma_eta_c = sqrt(par.var_eta_c);        
    
    % b. income variances
    par.var_xi  = max(-yy_lead1 - par.sigma_eta_y^2,0);
    par.var_psi = max(yy_perm,0);
        
        % standard deviations
        par.sigma_xi  = sqrt(par.var_xi);
        par.sigma_psi = sqrt(par.var_psi);
    
    % c. q and sigma_eps
    q       = -par.R*cy_lead1 -(par.R-1)*par.var_xi;
    q_limit = (sqrt(par.var_xi/par.var_psi+.25)-.5)*par.var_psi;
    
    if min([q par.var_xi par.var_psi]) <= 0
        par.sigma_eps = 0;
    elseif q >= q_limit
        par.sigma_eps = nan;
    else
        par.sigma_eps = sqrt( q*par.var_xi*(q+par.var_psi)/(par.var_xi*par.var_psi-q*(q+par.var_psi)) );
    end
    
end

end
end