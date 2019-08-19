clear;
clc;
close all;
delete('log_*')
funs.layout();

LOAD = 0;
do_analyze = 1;
do_profile = 1;
do_se = 1;
do_test = 1;

%% main

names = {'PT','pers','full'};

for i = 1:numel(names)
    
    name = names{i};
    estimate.run_and_analyze(name,'',{},[],do_analyze,do_profile,do_se,do_test,LOAD);

end