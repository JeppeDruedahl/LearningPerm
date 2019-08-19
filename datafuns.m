classdef datafuns
methods(Static)
    
function data = load(par)
    
    % 1. load
    import_age  = importdata(sprintf('psid/age%d.txt',par.group));
    import_C    = importdata(sprintf('psid/cons_res%d.txt',par.group));    
    import_Y    = importdata(sprintf('psid/income_res%d.txt',par.group));
    import_year = ones(size(import_age,1),1)*(1978:1:1992);
        
        % number of households
        N = size(import_age,1);    
        data.N = ceil(N/4)*4; % ensure modulo 4
    
    % 2. log
    import_logC = log(import_C);
    import_logY = log(import_Y);

    % 3. full age structure
    ages = par.age_min+(1:par.T);
        
    data.logC        = nan(data.N,par.T);
    data.logY        = nan(data.N,par.T);
    data.year        = nan(data.N,par.T);
    data.consecutive = nan(data.N,par.T);    
    
    for i = 1:N
       
       % a. unpack
       I = isnan(import_age(i,:)) == 0;
       age = import_age(i,I);
       logC = import_logC(i,I);
       logY = import_logY(i,I);    
       year = import_year(i,I);
       
       % b. indices
       J = max(ages == age',[],1);
       
       % c. pack
       data.logC(i,J) = logC;
       data.logY(i,J) = logY;
       data.year(i,J) = year;
              
    end
            
    % 4. analyze
    data = datafuns.analyze(data,par);
    
end
function data = analyze(data,par)
             
        % calculate min og max t
        ts = 1:par.T;
        I = sum(~isnan(data.logC) & ~isnan(data.logY)) > 0;
        data.min_t = min(ts(I));
        data.max_t = max(ts(I));
        
    % 1. moments
    if par.Nbootstraps > 0
        data_raw = data;
    else
        data_raw = [];
    end
    data = model.calc_moments(par,data,data.min_t,data.max_t);
               
    % 2. bootstrap
    data = datafuns.bootstrap_data(par,data,data_raw);    
        
    % 3. weight matrices
    if par.Nbootstraps > 0
        data.covmat = nancov(data.targets_vec_bs');
    end
    if par.do_full_output
        if par.Nbootstraps > 0
            data.covmat_bs = data.covmat;
        end
        data.covmat = datafuns.calc_covmat(par,data);
    end    
    data.W_full = data.covmat;
    data.W_diag = diag(diag(data.covmat));
    data.W_I = eye(numel(data.targets_vec))/1000; % numerically more stable with /1000
    
end
function data = select_func(par,data)

    if par.select_year > 0    
        data = datafuns.select_year_func(data,par.select_year);
    end
    if par.select_commault
        data = datafuns.select_commault_func(par,data);
    end
    
end
function data = select_year_func(par,data)

    I = data.year == par.select_year;
    vars = fieldnames(data);
    for j = 1:numel(vars)
       if contains(vars{j},'logY') || contains(vars{j},'logC')
            data.(vars{j})(~I) = nan;
       end
    end
        
end
function data = select_commault_func(par,data)
    
    I = ~isnan(data.DlogY) & ~isnan(data.DlogC) &...
        ~isnan(data.DlogY_lag1) & ~isnan(data.DlogY_lag2) &...
        ~isnan(data.DlogY_lead1) & ~isnan(data.DlogY_lead2);
    vars = fieldnames(data);
    for j = 1:numel(vars)
       if contains(vars{j},'logY') || contains(vars{j},'logC')
            data.(vars{j})(~I) = nan;
       end
    end
    
end
function data = bootstrap_data(par,data,data_raw)
        
        % not neeeded and par not output
        par.do_full_output = 0;
    
    % 1. field names
    moms_names = fieldnames(data.moms);
    moms_names = moms_names(~contains(moms_names,'_vec'));
    
    % 2. allocate
    for j = 1:numel(moms_names)
        data.moms.(sprintf('%s_bs',moms_names{j})) = nan(par.Nbootstraps,1);
    end
    data.mean_logC_bs = nan(size(data.logC,2),par.Nbootstraps);
    data.mean_logY_bs = nan(size(data.logY,2),par.Nbootstraps);
    data.targets_vec_bs = nan(numel(data.targets_vec),par.Nbootstraps);
    
        if par.Nbootstraps == 0
            return;
        end
    
    % 3. bootstrap
    for i = 1:par.Nbootstraps

        % a. dra
        bs_data = datafuns.bootstrap_sample(data_raw);

        % b. moments
        bs_data = model.calc_moments(par,bs_data,data.min_t,data.max_t,data.targets_I);       
        
        % c. save
        for j = 1:numel(moms_names)
            data.moms.(sprintf('%s_bs',moms_names{j}))(i) = bs_data.moms.(moms_names{j});
        end
        data.mean_logC_bs(:,i) = bs_data.mean_logC;
        data.mean_logY_bs(:,i) = bs_data.mean_logY;
        data.targets_vec_bs(:,i)  = bs_data.targets_vec;
        
    end

end
function bs_data = bootstrap_sample(data_raw)

    % a. ids
    [~,id]  = datasample(ones(data_raw.N,1),data_raw.N);
    id      = id';

    % b. pull
    vars = fieldnames(data_raw);
    for j = 1:numel(vars)
       if (contains(vars{j},'logY') || contains(vars{j},'logC')) && ...
               ~contains(vars{j},'mean')
            bs_data.(vars{j}) = data_raw.(vars{j})(id,:);
       end
    end
end
function covmat = calc_covmat(par,data)
    
    % a. build g matrix
    g = nan(numel(data.targets_vec),numel(data.logY(:)));
    j = 1;

    %fprintf('building gmat:\n')
    for i = 1:numel(par.targets)

        if strcmp(par.targets{i},'mean_logC') || strcmp(par.targets{i},'mean_logY')
            for t = 1:numel(data.(par.targets{i}))
                if t >= data.min_t && t <= min(par.max_max_t,data.max_t)
                if data.targets_I(i)
                    g(j,:) = data.(sprintf('%s_vec',par.targets{i}))(t,:);
                %fprintf('%12s, %2d: N = %6d: mean = %7.4f [%7.4f]\n',...
                %   par.targets{i},j,sum(~isnan(g(j,:))),nanmean(g(j,:)),data.(par.targets{i})(t));        
                j = j +1;
                end
                end
            end
        else
            if data.targets_I(i)            
                g(j,:) = data.moms.(sprintf('%s_vec',par.targets{i}));            
                %fprintf('%12s, %2d: N = %6d: mean = %7.4f [%7.4f]\n',...
                %    par.targets{i},j,sum(~isnan(g(j,:))),nanmean(g(j,:)),data.moms.(par.targets{i}));
                j = j + 1;
            end
        end    

    end
    
    % b. covmat
    Nmoms = numel(data.targets_vec);
    covmat = nan(Nmoms,Nmoms);
    for i = 1:Nmoms
    for j = 1:Nmoms
      
        % i. moments
        momi = data.targets_vec(i);
        momj = data.targets_vec(j);

        % ii. unit contributions
        gi = g(i,:);
        gj = g(j,:);

        % iii. not nan
        I = ~isnan(gi) & ~isnan(gj);
        
        % iv. number of obs
        Ni = sum(~isnan(gi));
        Nj = sum(~isnan(gj));

        % v. covmat
        covmat(i,j) = sum((gi(I)-momi).*(gj(I)-momj))/(Ni*Nj);

        %fprintf('(%2d,%2d): N = %6d: covmat = %12.8f\n',i,j,Ni*Nj,covmat(i,j));

    end
    end
    
    % c. remove nan
    I = isnan(covmat);
    covmat(I) = 0;

end


end
end