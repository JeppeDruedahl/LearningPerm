%% mex

intel = 1;
threads = 56;
compile_mex('mex_solve',threads,0,intel);
compile_mex('mex_solve',threads,1,intel);
compile_mex('mex_simulate',threads,0,intel);
compile_mex('mex_calc_logdiffs',threads,0,intel);
compile_mex('mex_calc_covs',threads,0,intel);

%% scripts

run_01_data
run_02_ceq
run_03_estimate
run_04_estimate_prefs
run_05_robustness
run_06_monte_carlo
run_07_tables
run_08_Commault