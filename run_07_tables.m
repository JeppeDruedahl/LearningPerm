clc;
clear;

%% 1. main

names = {'PT','pers','full','full_nocollege','full_college'};

pars = cell(numel(names),1);
for i = 1:numel(names)
    load(sprintf('data/%s',names{i}),'par');
    pars{i} = par;
end

postfix = '';
tabs.main_table(pars,postfix)

%% 2. ceq

names = {'ceq_measy2','ceq_group1_measy2','ceq_group2_measy2'};

pars = cell(numel(names),1);
for i = 1:numel(names)
    load(sprintf('data/%s',names{i}),'par');
    pars{i} = par;
end

postfix = '';
tabs.ceq(pars,postfix)

% low measurement error
names = {'ceq_measy1','ceq_group1_measy1','ceq_group2_measy1'};

pars = cell(numel(names),1);
for i = 1:numel(names)
    load(sprintf('data/%s',names{i}),'par');
    pars{i} = par;
end

postfix = '_tau_low';
tabs.ceq(pars,postfix)

% low measurement error
names = {'ceq_measy3','ceq_group1_measy3','ceq_group2_measy3'};

pars = cell(numel(names),1);
for i = 1:numel(names)
    load(sprintf('data/%s',names{i}),'par');
    pars{i} = par;
end

postfix = '_tau_high';
tabs.ceq(pars,postfix)

%% 3. robustness

namebases = {'PT','pers','full'};

for j = 1:numel(namebases)
    
namebase = namebases{j};    
    
specs = {{0,'rho_mid',{'rho'},2.0,'$\rho = 2.0$'},... 
         {0,'rho_high',{'rho'},4.0,'$\rho = 4.0$'},...  
         {0,'meas_y_frac_low',{'meas_y_frac'},0.0,'$\tau = 0.0$'},...  
         {0,'meas_y_frac_high',{'meas_y_frac'},0.50,'$\tau = 0.5$'},...  
         {0,'kappa_low',{'kappa'},0.25,'$\lambda = 0.25$'},...  
         {0,'kappa_high',{'kappa'},0.75,'$\lambda = 0.75$'},...  
         {0,'zeta_low',{'zeta'},{[0.2*ones(par.TR,1);zeros(par.T-par.TR,1)]},'$\zeta = 0.2$'},...       
         {0,'zeta_high',{'zeta'},{[0.4*ones(par.TR,1);zeros(par.T-par.TR,1)]},'$\zeta = 0.4$'},...                      
         {0,'W_diag',{'W_str'},{'W_diag'},'Diag. W'},...
         {0,'W_I',{'W_str'},{'W_I'},'Equal W'}};
     
pars = cell(numel(specs),1);
headers = cell(numel(specs),1);
for i = 1:numel(specs)
    
    name = sprintf('%s_%s',namebase,specs{i}{2}');
    load(sprintf('data/%s',name),'par');
    pars{i} = par;
    headers{i} = specs{i}{5};
    
end

postfix = sprintf('_%s',namebase);
tabs.robustness(pars,headers,postfix);

end

%% 4. robustness - borrowing constraint

namebases = {'PT','pers','full'};

for j = 1:numel(namebases)
    
namebase = namebases{j};    
    
specs = { {0,'zeta_low_alt',{'zeta'},{[0.10*ones(par.TR,1);zeros(par.T-par.TR,1)]},'$\zeta = 0.1$'},...
          {0,'zeta_low',{'zeta'},{[0.20*ones(par.TR,1);zeros(par.T-par.TR,1)]},'$\zeta = 0.2$'},...                          
          {0,'zeta_high_alt',{'zeta'},{[0.30*ones(par.TR,1);zeros(par.T-par.TR,1)]},'$\zeta = 0.3$'},...                              
          {0,'zeta_high',{'zeta'},{[0.40*ones(par.TR,1);zeros(par.T-par.TR,1)]},'$\zeta = 0.4$'}};
     
pars = cell(numel(specs)+1,1);
headers = cell(numel(specs)+1,1);
for i = 1:numel(specs)+1
    
    if i == 1
        load(sprintf('data/%s',namebase),'par');
        headers{i} = '';
    else
        name = sprintf('%s_%s',namebase,specs{i-1}{2}');
        headers{i} = specs{i-1}{5};
        load(sprintf('data/%s',name),'par');
    end
    pars{i} = par;
    
    
end

postfix = sprintf('_%s_borrowing_constraint',namebase);
tabs.robustness(pars,headers,postfix);

end