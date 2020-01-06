clear;
clc;
close all;
delete('log_*')
funs.layout();

LOAD = 0;
LOAD_PERS_ZERO = 0;
LOAD_PERS = 0;
LOAD_PERS_BETA = 0;
LOAD_PERS_PSI = 0;
LOAD_PERS_XI = 0;
LOAD_PERS_ALPHA = 0;

%% main

truepar = struct();
truepar.sigma_eps   = 0.10;
truepar.beta        = 0.96;
truepar.sigma_eta_c = 0.25;
truepar.sigma_eta_y = 0.10;
truepar.sigma_psi   = 0.15;
truepar.sigma_xi    = 0.15;
truepar.g0          = 0.06;
truepar.g1          = -0.05;
truepar.alpha       = 0.90;
truepar.omega       = 0.15;

modelname = 'full';
name = 'full';
N_MC = 5; % We only run 5 runs here because this code illustrates the setup but does NOT reproduce the results in the paper.

MC(modelname,name,N_MC,truepar,LOAD)
%close all;

%% pers

load('data/pers','par');
truepar = par;

modelname = 'pers';
name = 'pers_zero';
N_MC = 200;
MC(modelname,name,N_MC,truepar,LOAD_PERS_ZERO)
close all;

% with noise
truepar = par;
truepar.sigma_eps = 0.10;
name = 'pers';
MC(modelname,name,N_MC,truepar,LOAD_PERS)
close all;

% some noise + something more
truepar = par;
truepar.sigma_eps = 0.10;
truepar.beta = truepar.beta - 0.01;
name = 'pers_beta';
MC(modelname,name,N_MC,truepar,LOAD_PERS_BETA)
close all;

truepar = par;
truepar.sigma_eps = 0.10;
truepar.sigma_psi = sqrt(0.80)*truepar.sigma_psi;
name = 'pers_psi';
MC(modelname,name,N_MC,truepar,LOAD_PERS_PSI)
close all;

truepar = par;
truepar.sigma_eps = 0.10;
truepar.sigma_xi = sqrt(0.80)*truepar.sigma_xi;
name = 'pers_xi';
MC(modelname,name,N_MC,truepar,LOAD_PERS_XI)
close all;

truepar = par;
truepar.sigma_eps = 0.10;
truepar.alpha = 0.90;
name = 'pers_alpha';
MC(modelname,name,N_MC,truepar,LOAD_PERS_ALPHA)
close all;