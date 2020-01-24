# LearningPerm

Code for "[Can Consumers Distinguish Persistent from Transitory Income Shocks?](http://web.econ.ku.dk/druedahl/papers/2020_LearningPerm.pdf)", [Druedahl](http://web.econ.ku.dk/druedahl) and [Jørgensen](http://www.tjeconomics.com), 2020.

## Requirements

1. MATLAB 
2. C++ compiler

The code was run and tested with "MATLAB 2018b" and "Intel Parallel Studio XE 2018 for C++" on a 64-bit Windows machine.

**vectorclass:** The code relies on the **vectorclass** library developed by [Agnar Fog](https://www.agner.org/optimize/#vectorclass).

**Alternative C++ compiler:** You can alternatively install the "MinGW GCC" C++ compiler extension to MATLAB. Set `intel = 0` in `run_00_all.m`. If you use another MATLAB version change the path to `libgomp.a` in `compile_mex.m` accordingly. 

## ReadMe

Everything can be run from `run_00_all.m`. It calls all the `run_*.m` files in the correct order.

1. **Input:** PSID data in `psid/*.txt`.
2. **Output:** All figures and tables are saved in `fig_tabs/`.

## Reproduction

**Computation time:**  We have included a switch, `DO_SMALL`, in the main MATLAB file `run_00_all.m`. If its value is `0` the full Monte Carlo results with 200 runs will be produced. If its value is `1` (the default) only 5 Monte Carlo runs will be executed. Running the code with `DO_SMALL = 0` takes about 10 days and 12 hours using 56 threads on a computer with 2x Intel(R) Xeon(R) Gold 6154 3.00 GHz CPUs (18 cores, 36 logical processes each). Running the code with `DO_SMALL = 1` takes around 19 hours using the same computer specified above, but does not reproduce the Monte Carlo results in the paper.

**Mapping between folders and results in the paper:** All results are placed in the "figs_tab" folder which includes a set of sub-folder. Here is the list of figures and tables from the main text and the supplemental material:

1. Figure 1 and 2 is in the folder "full_MC"
1. Figure 3 is in the folder "full"
1. Figure 4 is in the folder "sigma_eps_full"
1. Figure 5 is in the folders "full_pref" (panels a and b) and "beta_full" (panels c and d)
1. Table 2 is "ceq.tex"
1. Table 3 is "main.tex"
1. Table 4 is "robustness_full.txt"
1. Figure C1 and C2 is in the folder "pers_zero_MC"
1. Figure C3 and C4 is in the folder "pers_MC"
1. Figure C5 is in the folder "pers_beta_MC" (panel a), "pers_psi_MC" (panel b), "pers_xi_MC" (panel c) and "pers_alpha_MC" (panel d)
1. Table C1 is "ceq_tau_low.tex"
1. Table C2 is "ceq_tau_high.tex"
1. Figure C6 is in the folder "ceq"
1. Figure C7 is in the folder "PT"
1. Figure C8 is in the folder "pers"
1. Figure C9 is in the folder "sigma_eps_PT"
1. Figure C10 is in the folder "sigma_eps_pers"
1. Table C3 is "robustness_PT.txt"
1. Table C4 is "robustness_pers.txt"
1. Figure C11 is in the folder "omega_Commault"
1. Figure C12 is in the folder "sigma_eps_PT" (panels a and b), "sigma_eps_pers" (panels c and d) and "sigma_eps_full" (panels e and f)

## PSID data

The used PSID data is found in `psid/*.txt`. 

To reproduce these txt-files run `psid/construct_data.do` varying the value of `scalar group` from 0 to 1 and 2. The code uses `psid/data3.dta` taken from the replication package for Blundell et. al. (2008). The `psid/construct_data.do` file is based on `mindist_AER.do` from this replication package.
