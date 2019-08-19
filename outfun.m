function stop = outfun(x,optimValues,state)
   
   global t0;

   stop = false;
   
   if ~strcmp(state,'init') && ~strcmp(state,'done')
        fprintf(' %3d: obj = %12.8f, sigma_eps = %12.8f',optimValues.iteration,optimValues.fval,x(1));
        if optimValues.iteration > 0
            fprintf(' (%4.1f secs)\n',toc(t0));
        else 
            fprintf('\n')
        end
   end
   t0 = tic;
   
   % check file to stop
%    fileID = fopen('stop_optimizer.txt','r');
%    stop = logical(cell2mat(textscan(fileID,'%d')));
%    fclose(fileID);
    
   % figure
%    switch state
%        case 'init'
%         
%            close all;
%            figure('name','fmincon')           
%            hold on           
% 
%            subplot(2,1,1)
%            grid on;
%            ylabel('fval')
%            
%            subplot(2,1,2)
%            grid on;           
%            ylabel('$\sigma_{\epsilon}$')
%            xlabel('iteration')
%            
%            shg;
%                       
%        case 'iter'               
%            
%            axes = get(gcf,'Children');
%            set(axes(1),'NextPlot','add')
%            set(axes(2),'NextPlot','add')
%           
%            plot(axes(2),optimValues.iteration,optimValues.fval,'o','Color','black','MarkerFaceColor','black');                  
%            plot(axes(1),optimValues.iteration,x(1),'o','Color','black','MarkerFaceColor','black');
%                    
%            xticks(axes(1),0:optimValues.iteration)
%            xticks(axes(2),0:optimValues.iteration)
%            
%            xlim(axes(1),[-0.1 optimValues.iteration+0.1])
%            xlim(axes(2),[-0.1 optimValues.iteration+0.1])           
%            
%            shg;
%            
%        case 'done'
%            
%            hold off
%        
%        otherwise
%            
%    end
%    
end