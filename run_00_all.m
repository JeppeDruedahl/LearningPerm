%% mex

intel = 1;
threads = 56;
compile_mex('mex_solve',threads,0,intel);
compile_mex('mex_solve',threads,1,intel);
compile_mex('mex_simulate',threads,0,intel);
compile_mex('mex_calc_logdiffs',threads,0,intel);
compile_mex('mex_calc_covs',threads,0,intel);

% time at start
file = fopen('timings_start.txt','w');
fprintf(file,'started at %s \n',datestr(datetime('now')));
fclose(file);
%% scripts

run_01_data
run_02_ceq
run_03_estimate
run_04_estimate_prefs
run_05_robustness

% Monte Carlo: choose baseline or small. _small is much faster (5 MC runs vs. 200) but does not reproduce results in the paper.
DO_SMALL = 0;
if DO_SMALL == 0
  run_06_monte_carlo
else
  run_06_monte_carlo_small
end

run_07_tables
run_08_Commault

% time at end
file = fopen('timings_end.txt','w');
fprintf(file,'ended at %s \n',datestr(datetime('now')));
fclose(file);
