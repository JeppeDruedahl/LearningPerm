clear;
clc;
close all;
delete('log_*')
funs.layout();

LOAD = 0;
do_analyze = 0;
do_profile = 0;
do_se = 0;
do_test = 0;

%% main

names = {'full'};

betas = 0.955:0.001:0.97;

if LOAD == 0
for h = 1:numel(names)
    
    name = names{h};
    load(sprintf('data/%s',name),'par');
    
for i = 1:numel(betas)
    
    % a. change paramters
    parnames = {'sigma_eps','beta','sigma_eta_c','sigma_xi','sigma_psi','g0','g1','alpha','omega','est_par'};
    vals = {par.sigma_eps,betas(i),par.sigma_eta_c,par.sigma_xi,par.sigma_psi,par.g0,par.g1,par.alpha,par.omega,{'sigma_eps'}};

    % b. postfix
    postfix = sprintf('beta%d',i);
    
    % c. run
    estimate.run_and_analyze(name,postfix,parnames,vals,do_analyze,do_profile,do_se,do_test,LOAD);
    
end
end
end

%% figures

for h = 1:numel(names)
       
% a. settings
name = names{h};

% b. load baseline
load(sprintf('data/%s',name),'par');
par_base = par;
par_base.figfolder = sprintf('%s_pref',name);

% c. load other
pars = cell(numel(betas),1);
alpha_vec = nan(numel(betas),1);
sigma_eps_vec = nan(numel(betas),1);
obj_vec = nan(numel(betas),1);
for i = 1:numel(betas)

	% a. load    	
	postfix = sprintf('beta%d',i);	
	load(sprintf('data/%s_%s',name,postfix),'par');

	% b. update
	pars{i} = par;
	pars{i}.displayname = sprintf('\beta = %5.3f',par.beta);

	% c. pull
    if par.sigma_eps < 2
        sigma_eps_vec(i) = par.sigma_eps;
        alpha_vec(i) = par.alpha;    
        obj_vec(i) = par.obj;
    else
        sigma_eps_vec(i) = nan;
        alpha_vec(i) = nan;    
        obj_vec(i) = nan;        
    end
end

% d. figure - obj
fig = figure('name','obj');
hold on;

	% i. base
	ax = plot(par_base.beta,par_base.obj,...
	    's','MarkerSize',10,'Linewidth', 1.5,'Color','black',...
	    'DisplayName','estimated $\beta$');   
                
	% ii. other
	ax = plot(betas,obj_vec,...
	    '-o','MarkerSize',5, 'Linewidth', 1.5,'Color',par.colors{1},...
	    'DisplayName','fixed $\beta$');   
	set(ax, 'MarkerFaceColor', get(ax, 'Color'));   

    % iii. details
    legend('show','Location','best')    
    xlabel('$\beta$')
    ylabel('objective function')
    grid on;

        ax = ancestor(ax,'axes');
        Yaxis = ax.YAxis;
        Yaxis.FontSize = 16;
        Xaxis = ax.XAxis;
        Xaxis.FontSize = 16;    
               
    % iv. lines
    min_val = min([obj_vec;par_base.obj])-10;
    max_val = max([obj_vec;par_base.obj])+10;
    ax = plot([par_base.beta par_base.beta],[min_val max_val],...
        '-','Color','black','Linewidth',1);
    funs.dont_display(ax);
    ylim([min_val max_val]);
            
    min_val = min([betas])-0.0025;
    max_val = max([betas])+0.0025;
    ax = plot([min_val max_val],[par_base.obj par_base.obj],...
        '-','Color','black','Linewidth',1);
    funs.dont_display(ax);
    xlim([min_val max_val]);
    
    % v. save           
    funs.printfig(par_base,fig);

% e. figure - sigma_eps
fig = figure('name','sigma_eps');
hold on;

	% i. base
	ax = plot(par_base.beta,par_base.sigma_eps,...
	    's','MarkerSize',10, 'Linewidth', 1.5,'Color','black',...
	    'DisplayName','estimated $\beta$');   

	% ii.. other
	ax = plot(betas,sigma_eps_vec,...
	    '-o','MarkerSize',5, 'Linewidth', 1.5,'Color',par.colors{1},...
	    'DisplayName','fixed $\beta$');   
	set(ax, 'MarkerFaceColor', get(ax, 'Color')); 

    % iii. details
    legend('show','Location','best')    
    xlabel('$\beta$')
    ylabel('$\sigma_{\epsilon}$')
    grid on;

        ax = ancestor(ax,'axes');
        Yaxis = ax.YAxis;
        Yaxis.FontSize = 16;
        Xaxis = ax.XAxis;
        Xaxis.FontSize = 16; 
               
    % iv. lines
    min_val = min([sigma_eps_vec;par_base.sigma_eps])-0.01;
    max_val = max([sigma_eps_vec;par_base.sigma_eps])+0.01;
    ax = plot([par_base.beta par_base.beta],[min_val max_val],...
        '-','Color','black','Linewidth',1);
    funs.dont_display(ax);
    ylim([min_val max_val]);
            
    min_val = min([betas])-0.0025;
    max_val = max([betas])+0.0025;
    ax = plot([min_val max_val],[par_base.sigma_eps par_base.sigma_eps],...
        '-','Color','black','Linewidth',1);
    funs.dont_display(ax);
    xlim([min_val max_val]);
    
    % v. save         
    funs.printfig(par_base,fig);

end

%% beta

for h = 1:numel(names)
       
% a. settings
name = names{h};

% b. load par
load(sprintf('data/%s',name),'par');

% c. setup analysis
vals = [0.96 0.9625 0.965 0.9675 0.97];   
vals_select = vals;

    select = max(vals == vals_select');
    varname = 'beta';
    postfix = name;
    latexvarname = '\beta';
    
    % data
    load('data/data','data')
    rng(data.rng_state);    
    [~,par.draws] = model.simulate(par,[]); 
        
model.analyze(par,data,varname,vals,select,postfix,latexvarname,-1);

end
