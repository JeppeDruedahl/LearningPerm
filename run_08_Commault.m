clear;
clc;
close all;
delete('log_*')

funs.layout();

LOAD = 0;

%% 1. analyze Commault

% a. load data
load('data/data','data');

% b. estimate model
load('data/PT','par');
par.omega = 0.15;
par.full_info = 1;
fprintf('%15s = %7.4f\n','rho',par.rho);
fprintf('%15s = %7.4f\n','beta',par.beta);
fprintf('%15s = %7.4f\n','rho',par.R);
fprintf('%15s = %7.4f\n','g0',par.g0);
fprintf('%15s = %7.4f\n','g1',par.g1);
fprintf('%15s = %7.4f\n','sigma_psi',par.sigma_psi);
fprintf('%15s = %7.4f\n','sigma_xi',par.sigma_xi);
fprintf('%15s = %7.4f\n','sigma_eta_y',par.sigma_eta_y);
fprintf('%15s = %7.4f\n','sigma_eta_c',par.sigma_eta_c);
fprintf('\n');

% b. update
parnames = {'Nphat','Nxihat','Niota','Na','Nm','Nsim'};
parvals = [50,31,5,100,300,200000];
par = funs.update_struct(par,parnames,parvals);

% c. analyze
vals = [0.00 0.05 0.10 0.15 0.20 0.25];

    select = vals == vals;
    varname = 'omega';
    postfix = 'Commault';
    latexvarname = '\omega';
    
model.analyze(par,data,varname,vals,select,postfix,latexvarname,LOAD);

% d. low beta 
par_alt = par;
par_alt.beta = 0.95;

postfix = 'Commault_lowbeta';
model.analyze(par_alt,data,varname,vals,select,postfix,latexvarname,LOAD);

% d. very low beta 
par_alt = par;
par_alt.beta = 0.94;

postfix = 'Commault_verylowbeta';
model.analyze(par_alt,data,varname,vals,select,postfix,latexvarname,LOAD);