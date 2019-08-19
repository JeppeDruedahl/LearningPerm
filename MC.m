function [] = MC(modelname,name,N_MC,truepar,LOAD)

    %% main

    % a. settings
    if modelname == 'full'
        est_par = {'sigma_eps','beta','sigma_eta_c','sigma_xi','sigma_psi','g0','g1','alpha','omega'};
        est_par_latex = {'\sigma_{\epsilon}','\beta','\sigma_c','\sigma_{\xi}','\sigma_{\psi}','g_0','g_1','\alpha','\omega'};
    elseif modelname == 'pers'
        est_par = {'sigma_eps','beta','sigma_eta_c','sigma_xi','sigma_psi','g0','g1','alpha'};
        est_par_latex = {'\sigma_{\epsilon}','\beta','\sigma_c','\sigma_{\xi}','\sigma_{\psi}','g_0','g_1','\alpha'};
    elseif modelname == 'PT'
        est_par = {'sigma_eps','beta','sigma_eta_c','sigma_xi','sigma_psi','g0','g1'};
        est_par_latex = {'\sigma_{\epsilon}','\beta','\sigma_c','\sigma_{\xi}','\sigma_{\psi}','g_0','g_1'};        
    end

    est_par_ceq = {'sigma_eps','sigma_eta_c','sigma_xi','sigma_psi'};
    est_par_ceq_latex = {'\sigma_{\epsilon}','\sigma_c','\sigma_{\xi}','\sigma_{\psi}'};

    % b. data
    load('data/data','data');
    rng(data.rng_state);

    % c. baseline and draws
    load(sprintf('data/%s',modelname),'par');
    if LOAD ~= 1  
        [~,draws] = model.simulate(par,[]);    
    end    
    simN = par.simN; % size of simulation in estimation
        
    % d. true
    true = par;
    true.simN        = 2000;

    true.sigma_eps   = truepar.sigma_eps;
    true.beta        = truepar.beta;
    true.sigma_eta_c = truepar.sigma_eta_c;
    true.sigma_eta_y = truepar.sigma_eta_y;
    true.sigma_psi   = truepar.sigma_psi;
    true.sigma_xi    = truepar.sigma_xi;
    true.g0          = truepar.g0;
    true.g1          = truepar.g1;
    true.alpha       = truepar.alpha;
    true.omega       = truepar.omega;
    true.meas_y_fra  = nan; % -> assume measurement error in income known

        % solve
        if LOAD ~= 1    
            [true, true_sol] = model.solve(true);
            save(sprintf('data/MC/%s_true',name),'true');        
        else        
            load(sprintf('data/MC/%s_true',name),'true');        
        end
        
    % e. allocate
    MC = struct();
    for i = 1:numel(est_par)
        MC.(est_par{i}) = nan(numel(est_par),1);
    end
    MC_ceq = struct();
    for i = 1:numel(est_par_ceq)
        MC_ceq.(est_par_ceq{i}) = nan(numel(est_par_ceq),1);
    end

    % f. loop
    for i_MC = 1:N_MC
        
        t0 = tic;
        fprintf('i_MC = %3d\n',i_MC);
        
        % a. data
        if LOAD ~= 1
            
            % i. simulate
            sim_MC = model.simulate(true,true_sol); % new draws used each time

            % ii. select
            vars = {'logC','logY'};

            for j = 1:numel(vars)

                data_MC.(vars{j}) = sim_MC.(vars{j});

                % age cutoffs
                data_MC.(vars{j})(:,1:data.min_t-1)   = nan;
                data_MC.(vars{j})(:,data.max_t+1:end) = nan;      

            end
            data_MC.N = true.simN;     
            
            % iii. analyze data
            true.do_full_output = 1;  
            true.Nbootstraps = true.Nbootstraps_ceq;
            data_MC = datafuns.analyze(data_MC,true);      
            true.do_full_output = 0;
            true.Nbootstraps = 0;
            
        end
        
        % b. setup from true and draws
        if LOAD ~= 1

            par = true;
            par.simN    = simN;
            par.est_par = est_par;
            par.draws   = draws; % same draws when estimating
                
        end
        
        % c. model: estimate or load
        filename = sprintf('data/MC/%s_MC_%d',name,i_MC);
        if LOAD == 0 || (LOAD ~= 1 && isfile(sprintf('%s.mat',filename)) == 0)
            [par,~,~,~] = estimate.run(par,data_MC);
            if par.sigma_eps == 2
                par.sigma_eps = 0;
                [par,~,~,~] = estimate.run(par,data_MC);                    
            end
            par.draws = [];
            save(filename,'par');    
        else
            load(filename,'par');
        end

        % d. ceq
        filename = sprintf('data/MC/%s_MC_ceq_%d',name,i_MC);
        if LOAD == 0 || (LOAD ~= 1 && isfile(sprintf('%s.mat',filename)) == 0)
            par_ceq = par;
            par_ceq.est_par = est_par_ceq;
            par_ceq = ceq.estimate(par,data_MC);        
            save(filename,'par_ceq');
        else
            load(filename,'par_ceq');
        end
        
        % e. save
        for i = 1:numel(est_par)
            MC.(est_par{i})(i_MC) = par.(est_par{i});
        end
        MC.obj(i_MC) = par.obj;   
        for i = 1:numel(est_par_ceq)
            MC_ceq.(est_par_ceq{i})(i_MC) = par_ceq.(est_par_ceq{i});
        end
        fprintf(' sigma_eps = %7.5f (time = %4.1f)\n',par.sigma_eps,toc(t0));
        
    end

    %% histograms - ceq

    true.figfolder = sprintf('%s_MC',name);

    for i = 1:numel(est_par_ceq)

        fig = figure('name',sprintf('ceq_%s',est_par_ceq{i}));
        
        ax = histogram(MC_ceq.(est_par_ceq{i}),15,...
                        'FaceColor',true.colors{1},...
                        'Normalization','probability','DisplayName','estimates');
        
        % lines
        v = axis; 
        line([true.(est_par_ceq{i}),true.(est_par_ceq{i})],v(3:4),'LineStyle','-',...
            'Color','black','LineWidth',2,'DisplayName','true');
        line([nanmean(MC_ceq.(est_par_ceq{i})) nanmean(MC_ceq.(est_par_ceq{i}))],v(3:4),'LineStyle','--',...
            'Color','black','LineWidth',2,'DisplayName','mean');
    
        str = sprintf('true: %5.3f\nmean: %5.3f\nstd.: %5.3f\nmedian: %5.3f',...
            true.(est_par_ceq{i}),nanmean(MC_ceq.(est_par_ceq{i})),...
            std(MC_ceq.(est_par_ceq{i})),...       
            prctile(MC_ceq.(est_par_ceq{i}),50));
        annotation('textbox',[0.65 0.6 0.3 0.3],'String',str,'FontSize',12,'FitBoxToText','on','Backgroundcolor','white');
        
        % details     
        if i == 1
            legend('show','Location','east')
        end
        xlabel(['$' est_par_ceq_latex{i} '$'],'FontSize',16)
        ylabel('density','FontSize',16)

            % axes
            ax = ancestor(ax,'axes');        
            Yaxis = ax.YAxis;
            Yaxis.FontSize = 16;
            Xaxis = ax.XAxis;
            Xaxis.FontSize = 16;

        ylim(v(3:4))
        v_x = axis; 
        xlim([v_x(1)-0.001 v_x(2)+0.001]);

        % e. save
        grid on;
        funs.printfig(true,fig);
        %close(fig);
        
    end

    %% histograms

    true.figfolder = sprintf('%s_MC',name);

    for i = 1:numel(est_par)

        fig = figure('name',sprintf('%s',est_par{i}));
        
        ax = histogram(MC.(est_par{i}),15,...
                        'FaceColor',true.colors{1},...
                        'Normalization','probability','DisplayName','estimates');        

        % lines
        v = axis; 
        line([true.(est_par{i}),true.(est_par{i})],v(3:4),'LineStyle','-',...
            'Color','black','LineWidth',1.5,'DisplayName','true');
        line([mean(MC.(est_par{i})) mean(MC.(est_par{i}))],v(3:4),'LineStyle','--',...
            'Color','black','LineWidth',1.5,'DisplayName','mean');
        
        str = sprintf('true: %5.3f\nmean: %5.3f\nstd.: %5.3f\nmedian: %5.3f',...
            true.(est_par{i}),mean(MC.(est_par{i})),...
            std(MC.(est_par{i})),...
            prctile(MC.(est_par{i}),50));
        annotation('textbox',[0.65 0.6 0.3 0.3],'String',str,'FontSize',12,'FitBoxToText','on','Backgroundcolor','white');
        
        % details     
        if i == 1
            legend('show','Location','east')
        end
        %xlabel(['$' est_par_latex{i} '$'],'FontSize',16)
        %ylabel('density','FontSize',16)
        
            % axes
            ax = ancestor(ax,'axes');        
            Yaxis = ax.YAxis;
            Yaxis.FontSize = 16;
            Xaxis = ax.XAxis;
            Xaxis.FontSize = 16;

        ylim(v(3:4))
        
        % e. save
        grid on;
        funs.printfig(true,fig);
        %close(fig);

end