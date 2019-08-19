classdef tabs
methods(Static)
    
function [] = main_table(pars,postfix)

    fileID = fopen(sprintf('figs_tabs\\main%s.tex',postfix),'w+');
    
    % 1. header
    fprintf(fileID,'\\begin{tabular}{clccccc} \\toprule \n'); 
    fprintf(fileID,' && \\multicolumn{3}{c}{Whole sample} & No college & College \\\\ \\cmidrule(lr){3-5} \\cmidrule(lr){6-6} \\cmidrule(lr){7-7} \n');
    fprintf(fileID,'  \\multicolumn{2}{l}{Parameter}& (1)& (2)& (3)& (4)& (5) \\\\ \\midrule \n');

    % 2. main rows
    rows = {{'sigma_eps','$\sigma_{\epsilon}$','Private signal (std.)',1},...
            {'beta','$\beta$','Discount factor',1},...
            {'sigma_eta_c','$\sigma_{c}$','Meas. error, cons. (std.)',1},...
            {'sigma_psi','$\sigma_{\psi}$','Persistent shock (std.)',1},...
            {'sigma_xi','$\sigma_{\xi}$','Transitory shock (std.)',1},...
            {'g0','$g_{0}$','Income growth, constant',1},...
            {'g1','$g_{1}$','Income growth, age',1},...
            {'alpha','$\alpha$','AR(1) component',1},...
            {'omega','$\omega$','MA(1) component',1}		
            };

    for i = 1:numel(rows)

        % a. estimates
        fprintf(fileID,' %s & %s',rows{i}{2},rows{i}{3});
        for j = 1:numel(pars)
            parname = rows{i}{1};
            if any(contains(pars{j}.est_par,parname))
                fprintf(fileID,' & $%4.3f$',pars{j}.(parname));
            else
                fprintf(fileID,' & $-$');		
            end
        end
        fprintf(fileID,'\\\\ \n');

        % b. se
        fprintf(fileID,' &');
        for j = 1:numel(pars)
            parname = rows{i}{1};
            if any(contains(pars{j}.est_par_se,parname))
                fprintf(fileID,' & $(%4.3f)$',pars{j}.se.(parname));
            else
                fprintf(fileID,' & $-$');		
            end
        end
        fprintf(fileID,'\\\\ \n');

    end
    
    % 3. objective and test
    fprintf(fileID,'\\midrule \\multicolumn{2}{l}{Objective}');
    for j = 1:numel(pars)
        fprintf(fileID,' & $%4.3f$',pars{j}.obj);
    end        
    fprintf(fileID,'\\\\ \n');
    
    fprintf(fileID,'\\multicolumn{2}{l}{p-value for $\\sigma_{\\epsilon}=0$}');
    for j = 1:numel(pars)
        if isfield(pars{j},'p_zero') && strcmp(pars{j}.W_str,'W_full')
            fprintf(fileID,' & $%4.3f$',pars{j}.p_zero);
        else
            fprintf(fileID,' & $-$');            
        end
    end
    
    
    % 4. footer    
    fprintf(fileID,'\\\\ \\bottomrule \n'); 
    fprintf(fileID,' \\end{tabular} \n');
    
    fclose(fileID);
    
end
function [] = ceq(pars,postfix)
    
    fileID = fopen(sprintf('figs_tabs\\ceq%s.tex',postfix),'w+');    
    
    % 1. header
    fprintf(fileID,'\\begin{tabular}{clccc} \\toprule \n');
    fprintf(fileID,'&& \\multicolumn{1}{c}{Whole sample} & No college & College \\\\   \\cmidrule(lr){3-3} \\cmidrule(lr){4-4} \\cmidrule(lr){5-5} \n');
    fprintf(fileID,' \\multicolumn{2}{l}{Parameter}& (1)& (2)& (3) \\\\ \\midrule \n');
    
    % 2. rows
    rows = {{'sigma_eps','$\sigma_{\epsilon}$','Private signal (std.)',1},...
            {'sigma_eta_c','$\sigma_{c}$','Meas. error, cons. (std.)',1},...
            {'sigma_psi','$\sigma_{\psi}$','Persistent shock (std.)',1},...
            {'sigma_xi','$\sigma_{\xi}$','Transitory shock (std.)',1},...
            };

    for i = 1:numel(rows)

        % a. estimates
        fprintf(fileID,' %s & %s',rows{i}{2},rows{i}{3});
        for j = 1:numel(pars)
            parname = rows{i}{1};
            fprintf(fileID,' & $%4.3f$',pars{j}.(parname));
        end
        fprintf(fileID,'\\\\ \n');

        % b. se
        fprintf(fileID,' &');
        for j = 1:numel(pars)
            parname = sprintf('%s_se',rows{i}{1});
            fprintf(fileID,' & $(%4.3f)$',pars{j}.(parname));            
        end
        fprintf(fileID,'\\\\ \n');

    end

    % 4. fotter
    fprintf(fileID,'\\bottomrule \n');     
    fprintf(fileID,'\\end{tabular} \n');

end
function [] = robustness(pars,headers,postfix)

    fileID = fopen(sprintf('figs_tabs\\robustness%s.tex',postfix),'w+');
    
    % 1. header
    fprintf(fileID,'\\begin{tabular}{cl');
    for i = 1:numel(pars)
        fprintf(fileID,'c');
    end
    fprintf(fileID,'} \\toprule \n'); 
    fprintf(fileID,'  & ');
    for i = 1:numel(pars)
        fprintf(fileID,'& %s',headers{i});
    end    
    fprintf(fileID,'\\\\ \n');    
    fprintf(fileID,'  \\multicolumn{2}{l}{Parameter}');
    for i = 1:numel(pars)
        fprintf(fileID,'& (%d)',i);
    end    
    fprintf(fileID,'\\\\ \\midrule \n');

    % 2. main rows
    rows = {{'sigma_eps','$\sigma_{\epsilon}$','Private signal (std.)',1},...
            {'beta','$\beta$','Discount factor',1},...
            {'sigma_eta_c','$\sigma_{c}$','Meas. error, cons. (std.)',1},...
            {'sigma_psi','$\sigma_{\psi}$','Persistent shock (std.)',1},...
            {'sigma_xi','$\sigma_{\xi}$','Transitory shock (std.)',1},...
            {'g0','$g_{0}$','Income growth, constant',1},...
            {'g1','$g_{1}$','Income growth, age',1},...
            {'alpha','$\alpha$','AR(1) component',1},...
            {'omega','$\omega$','MA(1) component',1}		
            };

    for i = 1:numel(rows)

        % a. estimates
        fprintf(fileID,' %s & %s',rows{i}{2},rows{i}{3});
        for j = 1:numel(pars)
            parname = rows{i}{1};
            if any(contains(pars{j}.est_par,parname))
                fprintf(fileID,' & $%4.3f$',pars{j}.(parname));
            else
                fprintf(fileID,' & $-$');		
            end
        end
        fprintf(fileID,'\\\\ \n');

        % b. se
        fprintf(fileID,' &');
        for j = 1:numel(pars)
            parname = rows{i}{1};
            if any(contains(pars{j}.est_par_se,parname))
                fprintf(fileID,' & $(%4.3f)$',pars{j}.se.(parname));
            else
                fprintf(fileID,' & $-$');		
            end
        end
        fprintf(fileID,'\\\\ \n');

    end
    
    % 3. objective and test
    fprintf(fileID,'\\midrule \\multicolumn{2}{l}{Objective}');
    for j = 1:numel(pars)
        fprintf(fileID,' & $%4.3f$',pars{j}.obj);
    end        
    fprintf(fileID,'\\\\ \n');
    
    fprintf(fileID,'\\multicolumn{2}{l}{p-value for $\\sigma_{\\epsilon}=0$}');
    for j = 1:numel(pars)
        if isfield(pars{j},'p_zero') && strcmp(pars{j}.W_str,'W_full')
            fprintf(fileID,' & $%4.3f$',pars{j}.p_zero);
        else
            fprintf(fileID,' & $-$');            
        end
    end
    
    % 4. footer    
    fprintf(fileID,'\\\\ \\bottomrule \n'); 
    fprintf(fileID,' \\end{tabular} \n');
    
    fclose(fileID);
    
end

end
end