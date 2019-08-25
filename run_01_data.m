clear;
clc;
close all;
delete('log_*')

%% 1. data

par = model.setup();

for group = [0 1 2] % 0: whole sample, 1: no college, 2: college %  1 2
for commault = [0]
for year = [0]
    
    % a. initialize random number generator
    rng(par.seed);
    
    % b. load
    par.group = group;
    par.select_commault = commault;
    par.select_year = year;    
    par.do_full_output = 1;
    par.Nbootstraps = 5000;
    data = datafuns.load(par);
   
        % save state of random number generator
        data.rng_state = rng; 
        
    % c. selection
    if year ~= 0
        fprintf('year = %4d\n',year);
    end
    if commault == 1
        fprintf('commault style\n')
    end
       
    % d. save
    filename = 'data/data';
    if group ~= 0
        filename = sprintf('%s_%d',filename,group);
    end    
    if commault == 1
        filename = sprintf('%s_commault',filename);
    end
    if year ~= 0
        filename = sprintf('%s_%d',filename,year);
    end
    save(filename,'data','-v7.3')
    
    % e. print
    if par.do_full_output
        SEs = diag(sqrt(diag(diag(data.covmat))));
    else
        SEs = nan(numel(data.targets_vec,1));
    end
    i = 1;
    for j = 1:size(data.targets_cell,1)
        
        % i. string
        target_str = data.targets_cell{j,2};

        % ii. pull
        if isfield(data,target_str)
            moms = data.(target_str);
            moms_bs = data.(sprintf('%s_bs',target_str));        
        else
            moms = data.moms.(target_str);
            moms_bs = data.moms.(sprintf('%s_bs',target_str));                    
        end
        
        % iii. print
        if strcmp(target_str,'mean_logC') || strcmp(target_str,'mean_logY')
            t_grid = data.min_t:min(par.max_max_t,data.max_t);
            for t = t_grid
                fprintf('%10s, t = %2d: %7.4f [se = %7.4f] [bs se = %7.4f]\n',target_str,t,...
                    moms(t),SEs(i),std(moms_bs(t,:)));
                i = i + 1;
            end
        else
            fprintf('%18s: %7.4f [se = %7.4f] [bs se = %7.4f]\n',target_str,moms,SEs(i),std(moms_bs));            
            i = i + 1;
        end
            
    end
    
    fprintf('%18s: %7d\n','observations',sum(isnan(data.logdiffs.DlogY(:)) == 0));
    fprintf('\n');
   
end
end
end