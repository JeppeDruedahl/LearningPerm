# LearningPerm

Code for "[Can Consumers Distinguish Persistent from Transitory Income Shocks?](http://web.econ.ku.dk/druedahl/papers/2018_LearningPerm.pdf)", [Druedahl](http://web.econ.ku.dk/druedahl) and [Jørgensen](http://www.tjeconomics.com), 2019.

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

**Computation time:** The computation time of all results in the paper was around 10 days and 12 hours, using 56 threads. Most of this time is spent generating the Monte Carlo section of the paper. 

We have included a switch, `DO_SMALL`, in the main Matlab file `run_00_all.m` where switching it to `1` (default is `0`) will reduce the number of MC runs from 200 (baseline, as in the paper) to 5. 
In turn, the baseline code with `DO_SMALL = 0` will reproduce the results in the paper but will take significant time to complete. Changing `DO_SMALL = 1` will run much faster (around 19 hours) but will not reproduce the MC results in the paper.  

## PSID data

The used PSID data is found in `psid/*.txt`. 

To reproduce these txt-files run `psid/construct_data.do` varying the value of `scalar group` from 0 to 1 and 2. The code uses `psid/data3.dta` taken from the replication package for Blundell et. al. (2008). The `psid/construct_data.do` file is based on `mindist_AER.do` from this replication package.
