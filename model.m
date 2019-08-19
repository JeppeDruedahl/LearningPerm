classdef model
methods(Static)

%%%%%%%%%%%%
% 1. setup %
%%%%%%%%%%%%

function par = setup()

    par = struct();
    par.seed = 35009210;
    
    % 1. demographics
    par.T    = 55;     % life-span
    par.TR   = 40;     % working-life-span
    
    par.tmin      = 0;   % first period to solve for    
    par.age_min   = 25;  % first age in simulations
    par.max_max_t = 35;  % do not use data after par.age_min+par.max_max_t for lcp's
    
    % 2. preferences
    par.rho  = 1.50;  % CRRA
    par.beta = 0.96;   % discount factor

    % 3. income
    par.g0        =  0.030; % constant
    par.g1        = -0.100;  % age
    par.g2        =  0.000;  % age^2
    
    par.sigma_xi  = 0.15; % transitory std.
    par.sigma_psi = 0.15; % permanent std.
    
    par.alpha     = 1.00;   % persistence
    par.omega     = 0.00;   % MA
    par.kappa     = 0.50;   % retirement
    
    % 4. information 
    par.sigma_eps  = 0.01;   % info std.       
    par.q0         = 0;      % initial information
    
    % 5. saving and borrowing
    par.R    = 1.03;                                        % interest rate factor
    par.zeta = [0.0*ones(par.TR,1);zeros(par.T-par.TR,1)];  % borrowing constraint
    par.A0   = 0.00;                                        % initial wealth in simulation

    % 6. grids         
    par.a_max = 15.0;   % max a value
    par.Na    = 100;    % number of grid points
    par.Niota = 6;      % number of quadrature points
        
        % additional grids
        par.Nm     = par.Na+1;
        par.Nphat  = 1;
        par.Nxihat = 1;
        
        par.phat_min  = -0.80;
        par.phat_max  = 2;
        par.xihat_min = -0.30;
        par.xihat_max = 0.30;        
        
    % 7. simulation
    par.simT = par.T;    % number of periods to simulate
    par.simN = 50000;    % number of households to simulate
    par.do_MPC_MPCP = 0; % calculate MPC and MPCP in simulations
    
    % 8. measurement error
    par.sigma_eta_y = 0.00;   % income
    par.sigma_eta_c = 0.27;   % consumption

    % 9. colors
    colors_hex = {'1f77b4','ff7f0e','2ca02c','d62728','9467bd',...
                  '8c564b','e377c2','7f7f7f','bcbd22','17becf'};
    par.colors = cell(numel(colors_hex),1);
    for i = 1:numel(colors_hex)
       par.colors{i} = funs.hex2rgb(colors_hex{i}); 
    end
    
    % 10. data and estimation
    par.group           = 0;
    par.select_commault = 0;
    par.select_year     = 0;
    par.do_full_output  = 0;    
    par.Nbootstraps     = 0;
    par.Nbootstraps_ceq = 5000;    
    par.Nbootstraps_est = 200;
    par.Nbootstraps_MC  = 200;    
    
    par.targets = {'mean_logC','mean_logY',...
                   'yy','yy_lead1','yy_lead2','yy_lead3','yy_lead4','yy_lead5',...
                   'cc','cc_lead1','cc_lead2','cc_lead3','cc_lead4','cc_lead5',...
                   'cy_lag2','cy_lag1','cy','cy_lead1','cy_lead2','cy_lead3','cy_lead4','cy_lead5'}; 
              
end
function par = create_grids(par)
            
    % 1. shocks

        % a. basic GaussHermite
        [iota_x, iota_w] = model.GaussHermite(1,0,par.Niota);
        
        iota_x_no = 0;
        iota_w_no = 1.0;
        
        % b. vectorize x    
        if par.Nxihat == 1

            [par.iota1_x,par.iota2_x,par.iota3_x] = ndgrid(iota_x,iota_x,iota_x);
            par.iota1_x = par.iota1_x(:);
            par.iota2_x = par.iota2_x(:);
            par.iota3_x = par.iota3_x(:);

            par.iota_x = [par.iota1_x par.iota2_x par.iota3_x]'; % used in MATLAB

            % c. vectorize weight     
            [par.iota1_w,par.iota2_w,par.iota3_w] = ndgrid(iota_w,iota_w,iota_w);
            par.iota1_w = par.iota1_w(:);
            par.iota2_w = par.iota2_w(:);
            par.iota3_w = par.iota3_w(:);

            par.iota_w = par.iota1_w .*par.iota2_w .* par.iota3_w; % used in MATLAB

        else
            
            if isfield(par,'full_info')
                if par.full_info == 1
                    [par.iota1_x,par.iota2_x,par.iota3_x,par.iota4_x,par.iota5_x] ...
                        = ndgrid(iota_x_no,iota_x_no,iota_x,iota_x,iota_x_no);
                else
                    [par.iota1_x,par.iota2_x,par.iota3_x,par.iota4_x,par.iota5_x] ...
                        = ndgrid(iota_x,iota_x,iota_x,iota_x,iota_x);                    
                end                
            else
                [par.iota1_x,par.iota2_x,par.iota3_x,par.iota4_x,par.iota5_x] ...
                    = ndgrid(iota_x,iota_x,iota_x,iota_x,iota_x);
            end
            
            par.iota1_x = par.iota1_x(:);
            par.iota2_x = par.iota2_x(:);
            par.iota3_x = par.iota3_x(:);
            par.iota4_x = par.iota4_x(:);
            par.iota5_x = par.iota5_x(:);

            % c. vectorize weight     
            if isfield(par,'full_info')
                if par.full_info == 1
                    [par.iota1_w,par.iota2_w,par.iota3_w,par.iota4_w,par.iota5_w] ...
                        = ndgrid(iota_w_no,iota_w_no,iota_w,iota_w,iota_w_no);
                else
                    [par.iota1_w,par.iota2_w,par.iota3_w,par.iota4_w,par.iota5_w] ...
                        = ndgrid(iota_w,iota_w,iota_w,iota_w,iota_w);                    
                end                
            else
                [par.iota1_w,par.iota2_w,par.iota3_w,par.iota4_w,par.iota5_w] ...
                    = ndgrid(iota_w,iota_w,iota_w,iota_w,iota_w);
            end            
            par.iota1_w = par.iota1_w(:);
            par.iota2_w = par.iota2_w(:);
            par.iota3_w = par.iota3_w(:);
            par.iota4_w = par.iota4_w(:);
            par.iota5_w = par.iota5_w(:);
            
            par.iota_w = par.iota1_w .*par.iota2_w .* par.iota3_w .* par.iota4_w .* par.iota5_w;
        
        end
        
        % d. count number of shock nodes
        par.Nshocks = numel(par.iota_w);

    % 3. grids
    par.grid_a = nan(par.Na,par.T);
    par.grid_m = nan(par.Nm,par.T);
    for t = 1:par.T
        par.grid_a(:,t) = linspace(-par.zeta(t)+1e-6,par.a_max,par.Na)';
        par.grid_m(:,t) = linspace(-par.zeta(t)+0,par.a_max,par.Nm)';
    end
    par.grid_phat  = linspace(par.phat_min,par.phat_max,par.Nphat);
    par.grid_xihat = linspace(par.xihat_min,par.xihat_max,par.Nxihat);
    
    % 4. calculate implied parameters
        
        % a. income growth
        agegrid      = (1:par.T)/100;         
        par.G        = exp(par.g0 + par.g1*(agegrid) + par.g2*(agegrid.^2));
        par.logG     = log(par.G);
        par.logkappa = log(par.kappa);
        
        % b. shocks
        par.var_psi = par.sigma_psi^2;
        par.mu_psi  = -0.5*par.var_psi;

        par.var_xi  = par.sigma_xi^2;   
        par.mu_xi   = -0.5*par.var_xi;

        % c. information       
        par.var_eps = par.sigma_eps^2;
        
        % b. measurement error
        par.var_eta_y = par.sigma_eta_y^2;
        par.var_eta_c = par.sigma_eta_c^2; 
    
    % 5. Kalman 
    par = model.Kalman(par);   
    
end
function [x, w] = GaussHermite(sigma,mu,n)

    i = 1:n-1;
    a = sqrt(i/2);
    CM = diag(a,1) + diag(a,-1);
    [V, L]   = eig(CM);
    [x, ind] = sort(diag(L));
    x = x*sqrt(2)*sigma + mu;
    V = V(:,ind)';
    w = (sqrt(pi) * V(:,1).^2)./sqrt(pi);
    
end
function par = Kalman(par)
          
    par.K = cell(par.T,1);
    par.q = nan(par.T,1);
    par.sigma_q_psi = nan(par.T,1);
        
        % omega > 0
        par.Q = cell(par.T,1);
        par.V = cell(par.T,1);
        par.D = cell(par.T,1);
        par.VsqrtD = cell(par.T,1);
        
    % 1. standard case
    if par.Nxihat == 1       
        
        for t = 1:par.T

            % a. var_eps == 0
            if par.var_eps == 0
                par.q(t) = 0;
                par.sigma_q_psi(t) = sqrt(par.var_psi);
                par.K{t} = [0 1];
                continue;
            end
        
            % b. var_eps < inf
            if par.var_eps < inf
            
                if t == 1
                    q_predict = par.alpha^2*par.q0 + par.var_psi;
                else 
                    q_predict = par.alpha^2*par.q(t-1,1) + par.var_psi;
                end     

                frac = q_predict / ( (par.var_xi+par.var_eps)*q_predict + par.var_xi*par.var_eps );

                par.K{t} = frac*[par.var_eps par.var_xi];           
                par.q(t) = (1-par.K{t}*[1;1])*q_predict;
            
            else

                if t == 1
                    q_predict = par.alpha^2*par.q0 + par.var_psi;
                else 
                    q_predict = par.alpha^2*par.q(t-1,1) + par.var_psi;
                end     

                par.K{t} = [q_predict / (q_predict + par.var_xi) 0];
                par.q(t) = (1-par.K{t}*[1;1])*q_predict;

            end % var_eps
            par.sigma_q_psi(t) = sqrt(par.alpha^2*par.q(t)+par.var_psi);
                    
        end  % t
    
    % 2. case with MA(1) term
    else
         
        assert(isinf(par.var_eps) == 0);
        
        % a. setup
        par.Fmat = [par.alpha,0,0;0,0,par.omega;0,0,0];
        Ft = transpose(par.Fmat);
        
        par.Wmat = [sqrt(par.var_psi),0,0;0,sqrt(par.var_xi),0;0,sqrt(par.var_xi),0];
        Wt = transpose(par.Wmat);
        
        par.Hmat = [1,1,0;1,0,0];
        Ht = transpose(par.Hmat);
        
        par.Rmat = [0,0;0,sqrt(par.var_eps)];
        Rt = transpose(par.Rmat);
        
        % b. Kalman recursion
        for t = 1:par.T
                    
            if t == 1
                Q_predict = par.Fmat*zeros(3,3)*Ft + par.Wmat*Wt;
            else
                Q_predict = par.Fmat*par.Q{t-1}*Ft + par.Wmat*Wt;
            end
            
            S = par.Hmat*Q_predict*Ht + par.Rmat*Rt;  
            par.K{t} = Q_predict*Ht/S;
            par.Q{t} = (eye(3,3)-par.K{t}*par.Hmat)*Q_predict;
            
            % eigenvalue decomposition
            [V,D] = eig(par.Q{t});
            [par.D{t},I] = sort(diag(D),'descend');
            par.V{t} = V(:, I);
            par.VsqrtD{t} = par.V{t}*sqrt(diag(par.D{t}));
             
            assert(abs(par.D{t}(3)) < 1e-12); % third eigenvalue is zero
            
        end % t
                
    end % omega        
    
end


%%%%%%%%%%%%
% 2. solve %
%%%%%%%%%%%%

function [par,sol]  = solve(par)

    par = model.create_grids(par); 
    if par.Nphat == 1
        sol = mex_solve_onedim(par);
    else
        sol = mex_solve(par);        
    end
    
end


%%%%%%%%%%%%%%%
% 3. simulate %
%%%%%%%%%%%%%%%

function [sim,draws] = simulate(par,sol,draws)
    
    % a. draw random shocks if not input
    if nargin < 3

        draws.q0    = normrnd(0,1,[par.simN,1]);        
        draws.xi    = normrnd(0,1,[par.simN,par.simT]);        
        draws.psi   = normrnd(0,1,[par.simN,par.simT]);
        draws.eps   = normrnd(0,1,[par.simN,par.simT]);
        draws.eta_y = normrnd(0,1,[par.simN,par.simT]);
        draws.eta_c = normrnd(0,1,[par.simN,par.simT]);
        
    end
    
        % return if no solution
        if numel(sol) == 0
            sim = struct;
            return 
        end
    
    % b. simulate
    sim = mex_simulate(par,sol,draws);

end


%%%%%%%%%%%%%%
% 4. moments %
%%%%%%%%%%%%%%

function sim = calc_moments(par,sim,min_t,max_t,targets_I)
    
    % a. log-differences
    %sim = model.calc_logdiffs(par,sim);    
    sim_demeaned.logC = sim.logC - nanmean(sim.logC,1);
    sim_demeaned.logY = sim.logY - nanmean(sim.logY,1);    
    
    sim.logdiffs = mex_calc_logdiffs(sim_demeaned);
        
    % b. select 
    sim = datafuns.select_func(par,sim);
    
    % c. covariances
    if isfield(par,'do_full_output')
        sim.moms = mex_calc_covs(sim.logdiffs,par.do_full_output);
    else
        sim.moms = mex_calc_covs(sim.logdiffs,0);        
    end
    
    % d. life-cycle profiles
    sim = model.calc_lcp(par,sim,min_t,max_t);
    
    % e. targets
    if nargin == 4
        [sim.targets_vec,sim.targets_I,sim.targets_cell] = estimate.targets(par,sim,min_t,max_t);
    elseif nargin == 5
        [sim.targets_vec,~,sim.targets_cell] = estimate.targets(par,sim,min_t,max_t,targets_I);
    else
        error('wrong number of inputs');
    end
    
end
function sim = calc_lcp(par,sim,min_t,max_t)
    
    % a. mean
    mean_logC = nanmean(sim.logC);
    mean_logY = nanmean(sim.logY);

    % b. allocate
    sim.mean_logC = nan(size(mean_logC));
    sim.mean_logY = nan(size(mean_logY));
    
    % c. demean
    t_grid = min_t:min(par.max_max_t,max_t);
    sim.mean_logC(t_grid) = mean_logC(t_grid) - nanmean(mean_logC(t_grid));
    sim.mean_logY(t_grid) = mean_logY(t_grid) - nanmean(mean_logY(t_grid));        
    
    % g. full output
    if isfield(par,'do_full_output')
    if par.do_full_output
        
        % i. demean unit data
        demean_logC = sim.logC - nanmean(mean_logC(t_grid));
        demean_logY = sim.logY - nanmean(mean_logY(t_grid));
        
        % ii. pack correctly
        mean_logC_vec = nan(numel(sim.mean_logC),size(sim.logC,1),size(sim.logC,2));
        mean_logY_vec = nan(numel(sim.mean_logY),size(sim.logY,1),size(sim.logY,2));               
        for t = 1:numel(sim.mean_logC)
            
            mean_logC_vec(t,:,t) = demean_logC(:,t);
            mean_logY_vec(t,:,t) = demean_logY(:,t);            
            
        end
        
        % iii. vectorize
        sim.mean_logC_vec = nan(numel(sim.mean_logC),numel(sim.logC(:)));
        sim.mean_logY_vec = nan(numel(sim.mean_logY),numel(sim.logY(:)));               
        for t = 1:numel(sim.mean_logC)
            sim.mean_logC_vec(t,:) = funs.vec((mean_logC_vec(t,:,:)));
            sim.mean_logY_vec(t,:) = funs.vec((mean_logY_vec(t,:,:)));    
        end
        
    end 
    end
    
end


%%%%%%%%%%
%%%%%%%%%%

function [] = fit_plots(par,sol,sim,data,postfix)
    
    % 1. income
    moms = {'yy','yy_lead1','yy_lead2','yy_lead3','yy_lead4','yy_lead5'};
    k_grid = 0:5;
    
    name = sprintf('covDyDy%s',postfix);
    x_label = 'k';    
    y_label = 'cov($\Delta y_t,\Delta y_{t+k}$)';

    model.basic_fit_plot(par,sim,data,moms,k_grid,name,x_label,y_label);
    
    % 2. consumption
    moms = {'cc','cc_lead1','cc_lead2','cc_lead3','cc_lead4','cc_lead5'};
    k_grid = 0:5;
    
    name = sprintf('covDcDc%s',postfix);
    x_label = 'k';    
    y_label = 'cov($\Delta c_t,\Delta c_{t+k}$)';

    model.basic_fit_plot(par,sim,data,moms,k_grid,name,x_label,y_label);

    % 3. mixed
    moms = {'cy_lag2','cy_lag1','cy','cy_lead1','cy_lead2','cy_lead3','cy_lead4','cy_lead5'};
    k_grid = -2:5;
    
    name = sprintf('covDcDy%s',postfix);
    x_label = 'k';    
    y_label = 'cov($\Delta c_t,\Delta y_{t+k}$)';

    model.basic_fit_plot(par,sim,data,moms,k_grid,name,x_label,y_label);
    
    % 4. life-cycle profiles
    if numel(sim) == 1
        sim.mean_Y = mean(exp(sim.logY),1);
        sim.mean_C = mean(exp(sim.logC),1);     
    else    
        for i = 1:numel(sim)
            sim{i}.mean_Y = mean(exp(sim{i}.logY),1);     
            sim{i}.mean_C = mean(exp(sim{i}.logC),1);     
        end
    end
    model.basic_fit_lcp_plot(par,sim,'mean_logC','$c_t$',data,'mean_logC',postfix);
    model.basic_fit_lcp_plot(par,sim,'mean_logY','$y_t$',data,'mean_logY',postfix);
    model.basic_fit_lcp_plot(par,sim,'mean_Y','$Y_t$',data,'mean_Y',postfix);
    model.basic_fit_lcp_plot(par,sim,'mean_C','$C_t$',data,'mean_C',postfix);
    
end
function [] = basic_fit_plot(par,sim,data,moms,k_grid,name,x_label,y_label)
    
    markers = {'o','s','d','+','x','^','v'};

    % 1. allocate
    sim_moms = nan(numel(moms),numel(par));
    data_moms = nan(numel(moms),1);
    data_moms_down = nan(numel(moms),1);
    data_moms_up = nan(numel(moms),1);    
    
    % 2. calculate
    for k = 1:numel(moms)
        
        % a. data
        data_moms(k) = data.moms.(moms{k});           
        data_SE = std(data.moms.(sprintf('%s_bs',moms{k})));
        
            data_moms_down(k) = data_moms(k) - 1.96*data_SE;
            data_moms_up(k) = data_moms(k) + 1.96*data_SE;
        
        % b. sim
        if numel(par) == 1
            sim_moms(k,1) = sim.moms.(moms{k});
        else
            for h = 1:numel(par)
                sim_moms(k,h) = sim{h}.moms.(moms{k});
            end       
        end
    end
    
    % 3. figure
    fig = figure('name',name);        
    hold on

    % a. data
    plot(k_grid,data_moms,'o','Color','black','DisplayName','data');

        ax = plot(k_grid,data_moms_down,'--','Color','black');
        funs.dont_display(ax);
        ax = plot(k_grid,data_moms_up,'--','Color','black');
        funs.dont_display(ax);

    % b. sim
    if numel(par) == 1
        ax = plot(k_grid,sim_moms,'-o','Color',par.colors{1},...
            'MarkerFaceColor',par.colors{1},'DisplayName','model');
    else
        for h = 1:numel(par)            
            ax = plot(k_grid,sim_moms(:,h),markers{h},'Color',par{h}.colors{h},...
                'MarkerFaceColor',par{h}.colors{h},'DisplayName',par{h}.displayname);
        end        
    end

    % c. x
    xlim([k_grid(1) k_grid(end)])
    set(gca,'XTick',k_grid)

        % axes
        ax = ancestor(ax,'axes');        
        Yaxis = ax.YAxis;
        Yaxis.FontSize = 16;
        Xaxis = ax.XAxis;
        Xaxis.FontSize = 16;    
        xlabel(x_label,'FontSize',16)
        ylabel(y_label,'FontSize',16)

    % d. legend
    legend('show','Location','best')

    % e. save
    grid on;
    if numel(par) == 1
        funs.printfig(par,fig);
    else
        funs.printfig(par{1},fig);        
    end  
    close(fig);
        
end
function [] = basic_fit_lcp_plot(par,sim,varname,latexname,data,datavarname,postfix)

    markers = {'o','s','d','+','x','^','v'};

    fig = figure('name',sprintf('%s%s',varname,postfix));
    hold on;
    
    % a. simulation
    if numel(par) == 1
        age_grid = par.age_min+(1:par.simT);
        ax = plot(age_grid,sim.(varname),...
            '-o','MarkerSize',5, 'Linewidth', 1.5,'Color',par.colors{1},...
            'MarkerFaceColor',par.colors{1},'DisplayName','model');       
    else
        age_grid = par{1}.age_min+(1:par{1}.simT);
        for h = 1:numel(par)
            linespec = sprintf('-%s',markers{h});
            ax = plot(par{h}.age_min+(1:par{h}.simT),sim{h}.(varname),...
                linespec,'MarkerSize',5, 'Linewidth', 1.5,'Color',par{h}.colors{h},...
                'MarkerFaceColor',par{h}.colors{h},'DisplayName',par{h}.displayname);       
        end
    end

    % b. data
    if isfield(data,datavarname)

        data_y = data.(datavarname);
        data_SE = nanstd(data.(sprintf('%s_bs',datavarname)),0,2)';

        plot(age_grid,data_y,'o','Color','black','DisplayName','data');        

        ax1 = plot(age_grid,data_y+1.96*data_SE,...
            '--','Color','black');
        funs.dont_display(ax1);

        ax2 = plot(age_grid,data_y-1.96*data_SE,...
            '--','Color','black');
        funs.dont_display(ax2);

    end

    % c. legend
    legend('show','Location','best')

        % axes
        ax = ancestor(ax,'axes');        
        Yaxis = ax.YAxis;
        Yaxis.FontSize = 16;
        Xaxis = ax.XAxis;
        Xaxis.FontSize = 16;    
        xlabel('age','FontSize',16)
        ylabel(latexname,'FontSize',16)

    % d. save
    grid on;
    if numel(par) == 1
        funs.printfig(par,fig);
    else
        funs.printfig(par{1},fig);        
    end    
    close(fig);
    
end


%%%%%%%%%%%%%%
% 6. analyze %
%%%%%%%%%%%%%%

function [] = analyze(par,data,varname,vals,select,postfix,latexvarname,LOAD)
    
    % 1. name
    if strcmp(postfix,'') == 0
        name = sprintf('%s_%s',varname,postfix);
    else
        name = varname;
    end
    fprintf('Running case %s\n',name);
    
        % clean up beforehand
        delete(sprintf('figs_tabs/%s/*',name));

    % 2. solve
    pars = cell(numel(vals),1);
    sols = cell(numel(vals),1);
    if LOAD ~= 1

        for i = 1:numel(vals)

            % a. setup and update
            if numel(par) == 1
                pars{i} = par;
                pars{i}.(varname) = vals(i);
                pars{i}.zeta = [pars{i}.zeta(1)*ones(pars{i}.TR,1);zeros(pars{i}.T-pars{i}.TR,1)];
                if isinf(vals(i)) || (strcmp(varname,'sigma_eps')) && vals(i) > 2
                    pars{i}.displayname = sprintf('$%s = \\infty$',latexvarname);
                elseif vals(i) == 0
                    pars{i}.displayname = sprintf('$%s = 0$',latexvarname);                
                else
                    pars{i}.displayname = sprintf('$%s = %5.3f$',latexvarname,vals(i));
                end                
            else
                pars{i} = par{i};
            end
            pars{i}.figfolder = name;            
            pars{i}.do_MPC_MPCP = 1;
            
            % b. solve
            t0 = tic;
            [pars{i},sols{i}] = model.solve(pars{i}); 
            
            fprintf(' sol:  %s = %5.3f (%3.1f secs)\n',varname,vals(i),toc(t0));
        
        end
    end
    
    % 3. save
    t0 = tic;
    filename = sprintf('data/%s.mat',name);
    if LOAD == 1
        load(filename,'vals','sols','pars');
        fprintf(' loaded in %3.1f secs\n',toc(t0));
    elseif LOAD == 0    
        save(filename,'vals','sols','pars','-v7.3');
        fprintf(' saved in %3.1f secs\n',toc(t0));
    end
      
    % 4. simulate and calculate moments
    sims = cell(numel(vals),1);
    for i = 1:numel(vals)
        
            % reset seed
            rng(data.rng_state)
        
        % a. simulate
        t0 = tic;
        sim = model.simulate(pars{i},sols{i});
        
        fprintf(' sim:  %s = %5.3f (%3.1f secs)\n',varname,vals(i),toc(t0));        
        
        % b. means
        lcp = struct;       
        lcp.mean_a    = mean(sim.A./exp(sim.p),1);      
        lcp.mean_MPC  = mean(sim.MPC,1);     
        lcp.mean_MPCP = mean(sim.MPCP,1); 
            
        % c. moments like in data
        vars = {'logC','logY','xi','psi','MPC','MPCP','A','p'};
        for j = 1:numel(vars)
            sims{i}.(vars{j}) = sim.(vars{j});
            sims{i}.(vars{j})(:,1:data.min_t-1)   = nan;
            sims{i}.(vars{j})(:,data.max_t+1:end) = nan;            
        end
        sims{i} = model.calc_moments(pars{i},sims{i},data.min_t,data.max_t,data.targets_I);
                    
            sims{i}.moms.MPC  = nanmean(sims{i}.MPC(:),1);     
            sims{i}.moms.MPCP = nanmean(sims{i}.MPCP(:),1); 
        
            % copy lcp
            lcp_vars = fieldnames(lcp);
            for j = 1:numel(lcp_vars)
                sims{i}.(lcp_vars{j}) = lcp.(lcp_vars{j});
            end
            
            % some checks
            omega = pars{i}.omega;
            var_psi = pars{i}.var_psi;
            var_xi = pars{i}.var_xi;
            var_eta_y = pars{i}.var_eta_y;
        
            fprintf('  check of the income process (assuming alpha = 1):\n');
            fprintf('   cov(dlogY,dlogY)  :    %8.5f (%8.5f in theory)\n',sims{i}.moms.yy,var_psi+(2+omega+2*omega^2)*var_xi+2*var_eta_y);
            fprintf('   cov(dlogY,dlogY+1):    %8.5f (%8.5f in theory)\n',sims{i}.moms.yy_lead1,(-1+2*omega-omega^2)*var_xi-var_eta_y);        
            fprintf('   cov(dlogY,dlogY+2):    %8.5f (%8.5f in theory) \n',sims{i}.moms.yy_lead2,-omega*var_xi);
            fprintf('   cov(dlogY,dlogY_perm): %8.5f (%8.5f in theory) \n',sims{i}.moms.yy_perm,-var_psi);
            
        % d. transitory transmission coefficient            
        trans_true = funs.cov(sims{i}.logdiffs.DlogC(:),sims{i}.xi(:))/var_xi;
        trans_BPP  = -sims{i}.moms.cy_lead1/((1-omega)*var_xi);
        trans_data = -data.moms.cy_lead1/((1-omega)*var_xi);
        trans_bias = omega/(1-omega)*funs.cov(funs.vec(sims{i}.logdiffs.DlogC(:,2:end)),funs.vec(sims{i}.xi(:,1:end-1)))/var_xi;
            
            % print
            fprintf('  transitory transmission coefficient:\n');
            fprintf('   MPC  : %7.4f\n',sims{i}.moms.MPC);
            fprintf('   true : %7.4f\n',trans_true);
            fprintf('   BPP  : %7.4f (data: %7.4f)\n',trans_BPP,trans_data);
            fprintf('   bias : %7.4f\n',trans_bias);

            % save
            sims{i}.moms.trans_true = trans_true;
            sims{i}.moms.trans_BPP  = trans_BPP;
            data.moms.trans_BPP     = trans_data;
            data.moms.trans_BPP_bs  = -data.moms.cy_lead1_bs/((1-omega)*var_xi);
            sims{i}.moms.trans_bias = trans_bias;
        
        % e. permanent transmission coefficient        
        perm_true = funs.cov(sims{i}.logdiffs.DlogC(:),sims{i}.psi(:))/var_psi;
        perm_BPP  = sims{i}.moms.cy_perm/var_psi;
        perm_data = data.moms.cy_perm/var_psi;        
        perm_bias = -funs.cov(funs.vec(sims{i}.logdiffs.DlogC(:,3:end)),funs.vec(sims{i}.xi(:,2:end-1))+omega*funs.vec(sims{i}.xi(:,1:end-2)))/var_psi;
            
            % print
            fprintf('  permaent transmission coefficient:\n');
            fprintf('   MPCP : %7.4f\n',sims{i}.moms.MPCP);
            fprintf('   true : %7.4f\n',perm_true');
            fprintf('   BPP  : %7.4f (data: %7.4f)\n',perm_BPP,perm_data);
            fprintf('   bias : %7.4f\n',perm_bias);

            % save
            sims{i}.moms.perm_true = perm_true;
            sims{i}.moms.perm_BPP  = perm_BPP;
            data.moms.perm_BPP     = perm_data;  
            data.moms.perm_BPP_bs  = data.moms.cy_perm_bs/var_psi;   
            sims{i}.moms.perm_bias = perm_bias;
            
    end
    
    % 5. figures
        
        % a. select     
        pars_fig = cell(sum(select),1);
        sims_fig = cell(sum(select),1);            
        sols_fig = cell(sum(select),1);        
        j = 1;
        for i = 1:numel(vals)
            if select(i)

                pars_fig{j} = pars{i};
                sims_fig{j} = sims{i};
                sols_fig{j} = sols{i};
                j = j + 1;
                
            end
        end            
        
        % b. moments and basic life-cycle        
        model.fit_plots(pars_fig,sols_fig,sims_fig,data,'');
        
        % c. more life-cycle
        vars = {{'mean_a','$A_t/P_t$',''},...
            {'mean_MPC','MPC',''},...
            {'mean_MPCP','MPCP',''}};
        
        for j = 1:numel(vars)

            lcp_varname = vars{j}{1};
            lcp_latexname = vars{j}{2};
            lcp_datavarname = vars{j}{3};

            model.basic_fit_lcp_plot(pars_fig,sims_fig,lcp_varname,lcp_latexname,data,lcp_datavarname,'');

        end
            
        % c. more moments    
        moments = {'yy','yy_lead1','yy_lead2','cc','cc_lead1','cy','cy_lead1',...
                   'MPC','trans_true','trans_BPP','trans_bias',...
                   'MPCP','trans_true','perm_BPP','perm_bias'};
           
        for j = 1:numel(moments)

            momname = moments{j};
            fig = figure('name',momname);

            % a. y-values
            y = nan(numel(vals),1);
            for i = 1:numel(vals)
                y(i) = sims{i}.moms.(momname);
            end
                    
            % b. sim
            if strcmp(varname,'sigma_eps')

                ax = plot(vals(2:end-1),y(2:end-1),...                
                    '-o','MarkerSize',5, 'Linewidth', 1.5,'Color',pars{1}.colors{1},...
                    'DisplayName','various $\sigma_\epsilon$');
                set(ax, 'MarkerFaceColor', get(ax, 'Color'));    
            
                % extremes
                hold on;
                K = numel(vals(1:end-1));
                plot(vals(1:end-1),y(1)*ones(K,1),...
                    '-v','Color',pars{1}.colors{1},...
                    'Linewidth',1,'DisplayName','$\sigma_\epsilon = 0$');
                plot(vals(1:end-1),y(end)*ones(K,1),...
                    '-^','Color',pars{1}.colors{1},...
                    'Linewidth',1,'DisplayName','$\sigma_\epsilon = \infty$');
                
            else

                ax = plot(vals,y,...
                    '-o','MarkerSize',5, 'Linewidth', 1.5,'Color',pars{1}.colors{1},...
                    'DisplayName','model');            
                hold on;
                set(ax, 'MarkerFaceColor', get(ax, 'Color'));    

            end
                       
            % c. data
            if strcmp(varname,'sigma_eps') && isfield(data.moms,momname)
                
                data_y    = data.moms.(momname);
                data_y_SE = std(data.moms.(sprintf('%s_bs',momname)));
                    
                    plot(vals(1:end-1),ones(numel(vals)-1,1)*data_y,...
                         '-','Linewidth', 1,'Color','black','DisplayName','data');                          

                    % confidence bounds
                    ax = plot(vals(1:end-1),ones(numel(vals)-1,1)*(data_y-1.96*data_y_SE),...
                        '--','Linewidth', 1,'Color','black');   
                    funs.dont_display(ax);
                    ax = plot(vals(1:end-1),ones(numel(vals)-1,1)*(data_y+1.96*data_y_SE),...
                        '--','Linewidth', 1,'Color','black');                   
                    funs.dont_display(ax);
                
                % legend
                legend('show','Location','best') 

            end
            
            % d. details
            grid on;

                ax = ancestor(ax,'axes');
                Yaxis = ax.YAxis;
                Yaxis.FontSize = 16;
                Xaxis = ax.XAxis;
                Xaxis.FontSize = 16;    
                xlabel(sprintf('$%s$',latexvarname),'FontSize',16)
                       
            % e. save           
            funs.printfig(pars{1},fig);
            close(fig);
        
        end
                
       % policy functions
       for i = 1:numel(pars_fig)
            model.policyfuncs(pars_fig{i},sols_fig{i},sprintf('_zeta%d',pars_fig{i}.zeta(1)*100))
        end

    % special figure 1 for Commault note
    if strcmp(varname,'omega') == 1
        
        nvals = numel(vals);
        
        fig = figure('name','Commault');
        hold on;
        
            % a. y-values
            momname = 'trans_true';
            y_trans_true = nan(nvals,1);
            for i = 1:nvals
                y_trans_true(i) = sims{i}.moms.(momname);
            end
            momname = 'trans_BPP';
            y_trans_BPP = nan(nvals,1);
            for i = 1:nvals
                y_trans_BPP(i) = sims{i}.moms.(momname);
            end
            momname = 'MPC';
            y_MPC = nan(nvals,1);
            for i = 1:nvals
                y_MPC(i) = sims{i}.moms.(momname);
            end
            
            % b. sim
            ax = plot(vals(1:nvals),y_MPC,...
                '-o','MarkerSize',5, 'Linewidth', 1.5,'Color',pars{1}.colors{1},...
                'DisplayName','MPC, model');            
            set(ax, 'MarkerFaceColor', get(ax, 'Color'));
            
            ax = plot(vals(1:nvals),y_trans_true,...
                '-o','MarkerSize',5, 'Linewidth', 1.5,'Color',pars{1}.colors{2},...
                'DisplayName','$\phi_{\xi} = $ cov$(\Delta c_{t},\xi_{t}) / \sigma_{\xi}^{2}$, model');            
            set(ax, 'MarkerFaceColor', get(ax, 'Color')); 
            
            ax = plot(vals(1:nvals),y_trans_BPP,...
                '-o','MarkerSize',5, 'Linewidth', 1.5,'Color',pars{1}.colors{3},...
                'DisplayName','$\hat{\phi}_{\xi}^{BPP}$, model');            
            set(ax, 'MarkerFaceColor', get(ax, 'Color')); 
            
            % c. data
            momname = 'trans_BPP';
            if isfield(data.moms,momname)
                
                data_y    = data.moms.(momname);
                data_y_SE = std(data.moms.(sprintf('%s_bs',momname)));
                
                plot(vals(1:nvals),ones(nvals,1)*data_y,...
                     '-','Linewidth', 1,'Color','black','DisplayName','$\hat{\phi}_{\xi}^{BPP}$, data');                          
                
            end
            
                % legend
                legend('show','Location','best') 
                
            % d. details
            grid on;

                ax = ancestor(ax,'axes');
                Yaxis = ax.YAxis;
                Yaxis.FontSize = 16;
                Xaxis = ax.XAxis;
                Xaxis.FontSize = 16;    
                xlabel(sprintf('$%s$',latexvarname),'FontSize',16)
                       
            % e. save           
            funs.printfig(pars{1},fig);
            close(fig);
            
    end
        
    if strcmp(varname,'omega') == 1
        
        nvals = numel(vals);
        
        fig = figure('name','Commault_bias');
        hold on;
        
            % a. y-values
            momname = 'trans_bias';
            y_trans_bias = nan(nvals,1);
            for i = 1:nvals
                y_trans_bias(i) = sims{i}.moms.(momname);
            end
     
            y_trans_bias_tot = nan(nvals,1);
            for i = 1:nvals
                y_trans_bias_tot(i) = sims{i}.moms.('trans_BPP')-sims{i}.moms.('trans_true');
            end
            
            % b. sim
            ax = plot(vals(1:nvals),y_trans_bias_tot,...
                '-o','MarkerSize',5, 'Linewidth', 1.5,'Color',pars{1}.colors{4},...
                'DisplayName','$\hat{\phi}_{\xi}^{BPP} - \phi_{\xi}$');            
            set(ax, 'MarkerFaceColor', get(ax, 'Color'));
            
            ax = plot(vals(1:nvals),y_trans_bias,...
                '-o','MarkerSize',5, 'Linewidth', 1.5,'Color',pars{1}.colors{5},...
                'DisplayName','$\omega/(1-\omega)$ cov$(\Delta c_{t},\xi_{t-1}) / \sigma_{\xi}^{2}$');            
            set(ax, 'MarkerFaceColor', get(ax, 'Color')); 
                        
                % legend
                legend('show','Location','best') 
                
            % d. details
            grid on;

                ax = ancestor(ax,'axes');
                Yaxis = ax.YAxis;
                Yaxis.FontSize = 16;
                Xaxis = ax.XAxis;
                Xaxis.FontSize = 16;    
                xlabel(sprintf('$%s$',latexvarname),'FontSize',16)
                       
            % e. save           
            funs.printfig(pars{1},fig);
            close(fig);
            
    end               
    
    if strcmp(varname,'sigma_eps') == 0
        return;
    end
    
        % e. transitory vs. MPC
        nvals = 7;
        fig = figure('name','trans_vs_MPC');
        hold on;
        
            % a. y-values
            momname = 'trans_BPP';
            y_trans = nan(nvals,1);
            for i = 1:nvals
                y_trans(i) = sims{i}.moms.(momname);
            end
            momname = 'MPC';
            y_MPC = nan(nvals,1);
            for i = 1:nvals
                y_MPC(i) = sims{i}.moms.(momname);
            end
            
            % b. sim
            ax = plot(vals(1:nvals),y_MPC,...
                '-s','MarkerSize',5, 'Linewidth', 1.5,'Color',pars{1}.colors{2},...
                'DisplayName','MPC, model');            
            set(ax, 'MarkerFaceColor', get(ax, 'Color'));
            
            ax = plot(vals(1:nvals),y_trans,...
                '-o','MarkerSize',5, 'Linewidth', 1.5,'Color',pars{1}.colors{1},...
                'DisplayName','$\hat{\phi}_{\xi}$, model');            
            set(ax, 'MarkerFaceColor', get(ax, 'Color'));    
            
            
            % c. data
            momname = 'trans_BPP';
            if isfield(data.moms,momname)
                
                data_y    = data.moms.(momname);
                data_y_SE = std(data.moms.(sprintf('%s_bs',momname)));
                
                plot(vals(1:nvals),ones(nvals,1)*data_y,...
                     '-','Linewidth', 1,'Color','black','DisplayName','$\hat{\phi}_{\xi}$, data');                          
                
                % confidence bounds
                ax = plot(vals(1:nvals),ones(nvals,1)*(data_y-1.96*data_y_SE),...
                    '--','Linewidth', 1,'Color','black');   
                funs.dont_display(ax);
                ax = plot(vals(1:nvals),ones(nvals,1)*(data_y+1.96*data_y_SE),...
                    '--','Linewidth', 1,'Color','black');                   
                funs.dont_display(ax);
                

            end
            
                % legend
                legend('show','Location','best') 
                
            % d. details
            grid on;

                ax = ancestor(ax,'axes');
                Yaxis = ax.YAxis;
                Yaxis.FontSize = 16;
                Xaxis = ax.XAxis;
                Xaxis.FontSize = 16;    
                xlabel(sprintf('$%s$',latexvarname),'FontSize',16)
                       
            % e. save           
            funs.printfig(pars{1},fig);
            close(fig);
    
        % f. permanent vs. MPCP
        fig = figure('name','perm_vs_MPCP');
        hold on;
        
            % a. y-values
            momname = 'perm_BPP';
            y_trans = nan(nvals,1);
            for i = 1:nvals
                y_trans(i) = sims{i}.moms.(momname);
            end
            momname = 'MPCP';
            y_MPC = nan(nvals,1);
            for i = 1:nvals
                y_MPC(i) = sims{i}.moms.(momname);
            end
            
            % b. sim
            ax = plot(vals(1:nvals),y_MPC,...
                '-s','MarkerSize',5, 'Linewidth', 1.5,'Color',pars{1}.colors{2},...
                'DisplayName','MPCP, model');            
            set(ax, 'MarkerFaceColor', get(ax, 'Color'));
            
            ax = plot(vals(1:nvals),y_trans,...
                '-o','MarkerSize',5, 'Linewidth', 1.5,'Color',pars{1}.colors{1},...
                'DisplayName','$\hat{\phi}_{\psi}$, model');
            set(ax, 'MarkerFaceColor', get(ax, 'Color'));    
            
            
            % c. data
            momname = 'perm_BPP';
            if isfield(data.moms,momname)
                
                data_y    = data.moms.(momname);
                data_y_SE = std(data.moms.(sprintf('%s_bs',momname)));
                
                plot(vals(1:nvals),ones(nvals,1)*data_y,...
                     '-','Linewidth', 1,'Color','black','DisplayName','$\hat{\phi}_{\psi}$, data');                          
                
                % confidence bounds
                ax = plot(vals(1:nvals),ones(nvals,1)*(data_y-1.96*data_y_SE),...
                    '--','Linewidth', 1,'Color','black');   
                funs.dont_display(ax);
                ax = plot(vals(1:nvals),ones(nvals,1)*(data_y+1.96*data_y_SE),...
                    '--','Linewidth', 1,'Color','black');                   
                funs.dont_display(ax);
                

            end
            
                % legend
                legend('show','Location','best') 
                
            % d. details
            grid on;

                ax = ancestor(ax,'axes');
                Yaxis = ax.YAxis;
                Yaxis.FontSize = 16;
                Xaxis = ax.XAxis;
                Xaxis.FontSize = 16;    
                xlabel(sprintf('$%s$',latexvarname),'FontSize',16)
                       
            % e. save           
            funs.printfig(pars{1},fig);
            close(fig);
            
end

function [] = policyfuncs(par,sol,postfix)
    
    sizem = size(sol.m{1});
    for t = []
        
        if sizem(2) == 1 
            i_phat_vec = [1];
        else
            i_phat_vec = [20,30,40];
        end
        for i_phat = i_phat_vec

            if sizem(2) == 1
                name = sprintf('cfunc_%d%s',t,postfix);
            else
                name = sprintf('cfunc_%d_p%d%s',t,i_phat,postfix);
            end
            
            fig = figure('name',name);
            ax = plot(sol.m{t}(:,i_phat),sol.c{t}(:,i_phat),...
                '-o','MarkerSize',3, 'Linewidth', 1.5,'Color',par.colors{1});            
            hold on;
            set(ax, 'MarkerFaceColor', get(ax, 'Color'));    
            grid on;

            ax = ancestor(ax,'axes');
            Yaxis = ax.YAxis;
            Yaxis.FontSize = 16;
            Xaxis = ax.XAxis;
            Xaxis.FontSize = 16;    
            xlabel('cash-on-hand, $m_t$','FontSize',16)
            ylabel('consumption, $c_t$','FontSize',16)
            xlim([-1 5])   
            ylim([0 3])

            funs.printfig(par,fig);
            close(fig);

       end
       
    if sizem(2) > 1
        
        [phat,m] = meshgrid(par.grid_phat,par.grid_m(:,t));
        Phat = exp(phat);
        M = m.*exp(phat);
        C = sol.c{t}.*exp(phat);

        fig = figure('name',sprintf('cfunc_%d%s',t,postfix));     
        ax = surf(M,Phat,C); 
        hold on;
        grid on;

        ax = ancestor(ax,'axes');
        Yaxis = ax.YAxis;
        Yaxis.FontSize = 16;
        Xaxis = ax.XAxis;
        Xaxis.FontSize = 16;    
        Zaxis = ax.XAxis;
        Zaxis.FontSize = 16;            
        xlabel('$m_t \cdot exp(\hat{p}_t)$','FontSize',16)
        ylabel('$exp(\hat{p}_t)$','FontSize',16)
        zlabel('$c_t \cdot exp(\hat{p}_t)$','FontSize',16)

        funs.printfig(par,fig);
        close(fig);

    end    
    end
    

end

end
end