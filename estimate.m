classdef estimate
methods(Static)

%%%%%%%%%%%%
% 1. setup %
%%%%%%%%%%%%

function par = setup()
        
    par = model.setup();
    
    % a. paramters to estimate
    par.est_par = {'sigma_eps','beta','sigma_eta_c','sigma_xi','sigma_psi','g0','g1','alpha','omega'};

    % b. measurement error like Pistaferri and Meghir (2004)
    par.meas_y_frac = 0.25;
        
    % c. weighting
    par.W_str       = 'W_full';
            
    % d. optimizer
    par.options_fmincon = optimoptions('fmincon','Display','none',...
                                       'Algorithm','sqp',...
                                       'MaxFunctionEvaluations',3000,...
                                       'OptimalityTolerance',1e-6,...
                                       'StepTolerance',1e-6,...
                                       'OutputFcn',@outfun);
    
end
function [lb,ub] = set_bounds(est_par)
    
    lb = nan(numel(est_par),1);
    ub = nan(numel(est_par),1);
    for i = 1:numel(est_par)
        if strcmp(est_par{i},'beta') 
            lb(i) = 0.80;
            ub(i) = 1.05;
        elseif strcmp(est_par{i},'sigma_xi') || strcmp(est_par{i},'sigma_psi') 
            lb(i) = 0.01;
            ub(i) = 0.25;
        elseif strcmp(est_par{i},'sigma_eta_c') 
            lb(i) = 0.00;
            ub(i) = 0.50;  
        elseif strcmp(est_par{i},'sigma_eps')   
            lb(i) = 0.00;
            ub(i) = 2.00; 
        elseif strcmp(est_par{i},'g0') 
            lb(i) = 0.00;
            ub(i) = 0.25;
        elseif strcmp(est_par{i},'g1') 
            lb(i) = -0.20;
            ub(i) = 0.10;
        elseif strcmp(est_par{i},'alpha') 
            lb(i) = 0.50;
            ub(i) = 1.00;  
        elseif strcmp(est_par{i},'omega') 
            lb(i) = -0.50;
            ub(i) = 0.50;            
        else
            error('unknown parameter');
        end
    end
     
end
function [targets_vec,targets_I,targets_cell] = targets(par,data,min_t,max_t,targets_I)
   
    % a. load targets
    targets_cell = cell(numel(par.targets),2);    
    for i = 1:numel(par.targets)
        target = par.targets{i};
        targets_cell{i,2} = target; 
        if strcmp(target,'mean_logC') || strcmp(target,'mean_logY')
            t_grid = min_t:min(par.max_max_t,max_t);
            targets_cell{i,1} = funs.vec(data.(target)(t_grid));
        else
            targets_cell{i,1} = data.moms.(target);
        end
    end
    
    % b. vectorize
    targets_vec = cell2mat(targets_cell(:,1));    
    if nargin == 4
        targets_I = ~isnan(targets_vec);
    end
    targets_vec = targets_vec(targets_I);
        
end


%%%%%%%%%%%%%%%%%
% 2. estimation %
%%%%%%%%%%%%%%%%%

function [] = run_and_analyze(namebase,postfix,parnames,vals,do_analyze,do_profile,do_se,do_test,LOAD)

    if strcmp(postfix,'')
        name = namebase;
    else
        name = sprintf('%s_%s',namebase,postfix);
    end
    fprintf('estimating: %12s\n',name);

    % a. setup
    par = estimate.setup();

        if strcmp(namebase,'PT')
        
            par.est_par = {'sigma_eps','beta','sigma_eta_c','sigma_xi','sigma_psi','g0','g1'};
            par.est_par_se = {'beta','sigma_eta_c','sigma_xi','sigma_psi','g0','g1'};              
        
        elseif strcmp(namebase,'pers')
            
            par.est_par = {'sigma_eps','beta','sigma_eta_c','sigma_xi','sigma_psi','g0','g1','alpha'}; 
            par.est_par_se = {'beta','sigma_eta_c','sigma_xi','sigma_psi','g0','g1','alpha'};
                
                % adjust initial values
                par.sigma_psi =  0.18;
                par.alpha     =  0.85;       
                par.g0        =  0.08;
                par.g1        = -0.04;
                
                % adjust grids
                par.Nphat   = 50;
                
        elseif strcmp(namebase,'full')
            
            par.est_par = {'sigma_eps','beta','sigma_eta_c','sigma_xi','sigma_psi','g0','g1','alpha','omega'};  
            par.est_par_se = {'beta','sigma_eta_c','sigma_xi','sigma_psi','g0','g1','alpha','omega'};  
            
                % adjust initial values
                par.alpha =  0.85;
                par.omega =  0.15;
                par.g0    =  0.10;
                par.g1    = -0.05;
                
                % adjust grid
                par.Nphat   = 50;
                par.Nxihat  = 10;
                par.Niota   = 4;
                par.Na      = 80;
                par.Nm      = 150;
        
        else
            error('unknown model');
        end
            
    % b. update
    par = funs.update_struct(par,parnames,vals);
    par.zeta = [par.zeta(1)*ones(par.TR,1);zeros(par.T-par.TR,1)];
    
    % c. load data 
    filename = 'data/data';        
    if par.group ~= 0
        filename = sprintf('%s_%d',filename,par.group);            
    end
    load(filename,'data');          

    % d. estimate
    if LOAD == 0
       
        % i. draw
        rng(data.rng_state);    
        [~,par.draws] = model.simulate(par,[]); 
        
        % ii. estimate
        
            if do_profile
                profile clear;
                profile on;
            end
                
        if numel(par.est_par) > 0
            [par,sol,~,sim_data] = estimate.run(par,data);
        else
            [obj,par,sol,~,sim_data] = estimate.obj_fun(par,data,[]);
            par.obj = obj;
        end
        
            if do_profile
                profile off;
                profsave(profile('info'),sprintf('profiling/profile_%s',name))       
            end        
            
        
    end
    
    filename = sprintf('data/%s',name);
    if LOAD == 1
        load(filename,'par')
        rng(data.rng_state);    
        [~,par.draws] = model.simulate(par,[]);        
        [~,~,sol,~,sim_data] = estimate.obj_fun(par,data,[]);     
    end         
    
        % iii. standard errors   
        if (LOAD == 0 && do_se) || do_se == 2
            par = estimate.standard_errors(par,data);    
        end        

        % iv. test
        if (LOAD == 0 && do_test) || do_test == 2

            fprintf('test for sigma_eps = 0\n');

            par_h0 = par;
            par_h0.draws = par.draws;
            par_h0.sigma_eps = 0;
            par_h0.est_par = par_h0.est_par(2:end);
            [par_h0,~,~,~] = estimate.run(par_h0,data);

                qlr = par_h0.obj-par.obj;
                par.p_zero = 1-(0.5+0.5*chi2cdf(qlr,1));
                par.obj_zero = par_h0.obj;
                
            fprintf('\n')    

        end
    
    % e. save 
    if LOAD == 0 || do_se == 2 || do_test == 2
        draws = par.draws;
        par.draws = [];        
        save(filename,'par')
        par.draws = draws;
    end
    
    % f. parameters
    for j = 1:numel(par.est_par)
       fprintf('%12s: %7.4f\n',par.est_par{j},par.(par.est_par{j}));     
    end    
    fprintf('%12s: %7.4f\n','obj',par.obj);
    
    % g. fit figures
    par.figfolder = name;
    model.fit_plots(par,sol,sim_data,data,'')
    
    % h. analysis figures
    if do_analyze
        
        if par.Nxihat == 1
            vals = [0.00 0.01 0.025 0.05 0.10 0.25 0.50 0.75 1.0 1.5 2.0 inf];   
            vals_select = [0.00 0.05 0.10 0.25 inf];
        else
            vals = [0.00 0.01 0.025 0.05 0.10 0.25 0.50 0.75 1.0 1.5 2.0 5.0];   
            vals_select = [0.00 0.05 0.10 0.25]; 
        end

            select = max(vals == vals_select');
            varname = 'sigma_eps';
            postfix = name;
            latexvarname = '\sigma_{\epsilon}';

        model.analyze(par,data,varname,vals,select,postfix,latexvarname,-1);

    end
    fprintf('\n')        
    
    model.policyfuncs(par,sol,'')

end
function [par,sol,sim,sim_data] = run(par,data)
    
        fileID = fopen('stop_optimizer.txt','w+');
        fprintf(fileID,'0\n');
        fclose(fileID);

    % 1. grids
    assert((par.alpha == 1 && par.Nxihat == 1) || par.Nphat > 1);
    assert(par.omega == 0 || par.Nxihat > 1);

    % 2. initial
    theta0 = nan(numel(par.est_par),1);
    for p = 1:numel(par.est_par)
        theta0(p) = par.(par.est_par{p});
    end
    [lb,ub] = estimate.set_bounds(par.est_par);
        
    % 3. estimate
    obj_fun = @(theta) estimate.obj_fun(par,data,theta);
    [theta,par.obj,~] = fmincon(obj_fun,theta0,[],[],[],[],lb,ub,[],par.options_fmincon);
    %[theta,par.obj,~] = fminsearch(obj_fun,theta0,par.options_fminsearch);
    
    % 4. details
    [~,par,sol,sim,sim_data] = estimate.obj_fun(par,data,theta);
    
end
function [obj,par,sol,sim,sim_data,obj_diff_vec] = obj_fun(par,data,theta)

    % 1. update
    if numel(theta) > 0
        par = funs.update_struct(par,par.est_par,theta);
    end
    
    % 2. measurement
    if ~isnan(par.meas_y_frac)
        par.sigma_eta_y = sqrt(0.5*par.meas_y_frac*data.moms.yy);   
    end
    
    % 3. solve
    [par, sol] = model.solve(par);

    % 4. simulate
    sim = model.simulate(par,sol,par.draws);
    
    % 5. demeand and nanify like data
    vars = {'logC','logY','A',};
    for j = 1:numel(vars)
        sim_data.(vars{j}) = sim.(vars{j});
        sim_data.(vars{j})(:,1:data.min_t-1)   = nan;
        sim_data.(vars{j})(:,data.max_t+1:end) = nan;            
    end
        
    % 6. moments 
    sim_data = model.calc_moments(par,sim_data,data.min_t,data.max_t,data.targets_I);  

    % 7. objective function
    obj_diff_vec = sim_data.targets_vec - data.targets_vec;
    obj = obj_diff_vec'*(data.(par.W_str)\obj_diff_vec);
    
end
function par = standard_errors(par,data)
        
        % replace est_par
        est_par = par.est_par;
        par.est_par = par.est_par_se;
 
    % 1. setup
    num_mom     = numel(data.targets_vec);
    num_est_par = numel(par.est_par_se);

    % 2. theta
    theta = nan(num_est_par,1);
    for i = 1:num_est_par
        theta(i) = par.(par.est_par_se{i});
    end
    
    % 3. gradients
    h    = 1.0e-5;
    grad = nan(num_mom,num_est_par);
        
        [~,~,~,~,~,valbase] = estimate.obj_fun(par,data,theta);
        
        for i = 1:num_est_par
            
            % a. theta
            theta_now    = theta;                   
            theta_now(i) = theta(i)+h;
            
            % b. objective
            [~,~,~,~,~,forward] = estimate.obj_fun(par,data,theta_now);
            
            % c. gradient
            grad(:,i)    = -(forward-valbase)./(h);
            
        end

    % 2. matrices
    W = data.(par.W_str);    
    S = data.covmat; % data.Nmat(i,j) = N(i)*N(j); data.cov./(data.Nmat)
    
        GW   = grad'/W; % grad'*inv(W)
        GWG  = grad'*(W\grad); %grad'*inv(W)*grad = GW*grad
    
    % 3. standard errors
    Avar    = (GWG\(GW*S*GW'))/GWG;
    SE      = sqrt(diag(Avar));
    
    % 4. save
    for i = 1:num_est_par
        par.se.(par.est_par_se{i}) = SE(i);
    end
    
        % restpre est_par
        par.est_par = est_par;
        
end

end
end