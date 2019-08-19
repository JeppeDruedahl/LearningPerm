clear;
clc;
close all;
delete('log_*')

funs.layout();

%% main

est_par = {'sigma_eps','sigma_eta_c','sigma_xi','sigma_psi'};
est_par_latex = {'\sigma_{\epsilon}','\sigma_c','\sigma_{\xi}','\sigma_{\psi}'};
meas_y_frac = [0.0 0.25 0.50];
pars = cell(numel(meas_y_frac),1);

for group = [0 1 2]
    
% a. load data
if group == 0
    load('data/data','data');
else
    load(sprintf('data/data_%d',group),'data');    
end

% b. estimate
for i = 1:numel(meas_y_frac)
    
    % a. setup
    par = estimate.setup();
    
        par.est_par     = est_par;
        par.meas_y_frac = meas_y_frac(i);
    
    % b. estimate
    par = ceq.estimate(par,data);
    
    % c. save
    if group == 0
        filename = sprintf('data/ceq_measy%d',i);
    else
        filename = sprintf('data/ceq_group%d_measy%d',group,i);        
    end    
    save(filename,'par');
    
    if group == 0
        pars{i} = par;
    end
    
end
end


%% histogram

pars{1}.figfolder = 'ceq';

for j = 1:numel(est_par)

    fig = figure('name',sprintf('%s',est_par{j}));
    hold on;
    
    for i = 1:numel(meas_y_frac)    
        
        vals = pars{i}.(sprintf('%s_bs',est_par{j}));        
        I_nan = isnan(vals);
        I_zero = vals == 0; 
        
        if strcmp(est_par{j},'sigma_eps')
            displayname = sprintf('$\\tau = %3.2f$: Pr$[\\sigma_{\\epsilon}=0] = %3.2f$',...
                meas_y_frac(i),mean(I_zero));
        else
            displayname = sprintf('$\\tau = %3.2f$',meas_y_frac(i));
        end
 
        ax = histogram(vals,'BinWidth',0.002,...
            'FaceColor',pars{i}.colors{i},...
            'Normalization','probability','DisplayName',displayname);    
        
    end
    if strcmp(est_par{j},'sigma_eps')
        xlim([0 0.15])
        ylim([0 0.05])
    end 
    
    % lines
    v = axis;     
    for i = 1:numel(meas_y_frac)
        ax = line([pars{i}.(est_par{j}),pars{i}.(est_par{j})],v(3:4),'LineStyle','--',...
            'Color',pars{i}.colors{i},'LineWidth',1.5);
        funs.dont_display(ax);
    end
    
    % details  
    if strcmp(est_par{j},'sigma_eps')    
        legend('show','Location','northeast')
    else
        legend('show','Location','best')
    end
    xlabel(['$' est_par_latex{j} '$'],'FontSize',16)
    ylabel('density','FontSize',16)

        % axes
        ax = ancestor(ax,'axes');        
        Yaxis = ax.YAxis;
        Yaxis.FontSize = 16;
        Xaxis = ax.XAxis;
        Xaxis.FontSize = 16;
       
    ylim(v(3:4))
    
    % e. save
    grid on;
    funs.printfig(pars{1},fig);
    close(fig);
    
end