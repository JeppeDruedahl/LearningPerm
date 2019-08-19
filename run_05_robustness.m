clear;
clc;
close all;
delete('log_*')

funs.layout();

do_analyze = 0;
do_profile = 0;
do_se = 1;
do_test = 1;

%% main

names = {'PT','pers','full'};

par = estimate.setup();
specs = {{0,'nocollege',{'group'},1},...
         {0,'college',{'group'},2},... 
         {0,'rho_mid',{'rho'},2.0,'$\rho = 2.0$'},... 
         {0,'rho_high',{'rho'},4.0,'$\rho = 4.0$'},...  
         {0,'meas_y_frac_low',{'meas_y_frac'},0.0,'$\tau = 0.0$'},...  
         {0,'meas_y_frac_high',{'meas_y_frac'},0.50,'$\tau = 0.5$'},...  
         {0,'kappa_low',{'kappa'},0.25,'$\lambda = 0.25$'},...  
         {0,'kappa_high',{'kappa'},0.75,'$\lambda = 0.75$'},...  
         {0,'zeta_low_alt',{'zeta'},{[0.10*ones(par.TR,1);zeros(par.T-par.TR,1)]},'$\zeta = 0.1$'},...
         {0,'zeta_low',{'zeta'},{[0.20*ones(par.TR,1);zeros(par.T-par.TR,1)]},'$\zeta = 0.2$'},...                          
         {0,'zeta_high_alt',{'zeta'},{[0.30*ones(par.TR,1);zeros(par.T-par.TR,1)]},'$\zeta = 0.3$'},...                              
         {0,'zeta_high',{'zeta'},{[0.40*ones(par.TR,1);zeros(par.T-par.TR,1)]},'$\zeta = 0.4$'},...                              
         {0,'W_diag',{'W_str'},{'W_diag'},'Diag. W'},...
         {0,'W_I',{'W_str'},{'W_I'},'Equal W'}};
     
for i = 1:numel(names)
for j = 1:numel(specs)
    
    % a. name
    name = names{i};
    
    % b. load
    LOAD = specs{j}{1};
    if LOAD == -1
        continue
    end
    
    % c. settings
	postfix  = specs{j}{2};
	parnames = specs{j}{3};
	parvals  = specs{j}{4};	

    % d. run
	estimate.run_and_analyze(name,postfix,parnames,parvals,do_analyze,do_profile,do_se,do_test,LOAD);
    
end
end