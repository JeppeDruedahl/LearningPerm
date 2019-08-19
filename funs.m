classdef funs
methods(Static)

%%%%%%%%%%%%%
% 1. output %
%%%%%%%%%%%%%

function [] = layout()

    % set layout parameters
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultTextInterpreter','latex');
    set(groot, 'defaultAxesFontSize', 12); 
    warning('off','all')

end  
function RGB = hex2rgb(HEX)

	if ~iscellstr(HEX)
		HEX = cellstr(HEX);
	end
	HEX = upper(HEX);
	HEX = double(char(HEX));
	tf = HEX > 64;
	HEX(~tf) = HEX(~tf) - '0';
	HEX(tf) = HEX(tf) - 'A' + 10;

	HEX = reshape(HEX.',2,[]);
	HEX(1,:) = 16*HEX(1,:);
	RGB = sum(HEX,1);
	RGB  = reshape(RGB,3,[]);
	RGB = RGB.';
    RGB = RGB/255;

end
function [] = printfig(par,figin)

    fig = figure(figin);
    fig.PaperUnits = 'centimeters';   
    fig.PaperPositionMode = 'manual';
    fig.PaperPosition = [0 0 16 12];
    fig.PaperSize = [16 12];

    title('');
    
    folder = ['figs_tabs\' par.figfolder];
    if exist(folder,'dir') == 0
        mkdir(folder);
    end
    filename = [folder '\' get(figin,'name') ''];
    print('-dpdf',['' filename '.pdf']);
        
end
function dont_display(h)
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end  

%%%%%%%%%%%
% 2. misc %
%%%%%%%%%%%

function y = vec(x)
    y = x(:);
end  
function x = nonlinspace(lo,hi,n,phi)
    % recursively constructs an unequally spaced grid.
    % phi > 1 -> more mass at the lower end of the grid.
    % lo can be a vector (x then becomes a matrix).

    x      = NaN(n,length(lo));
    x(1,:) = lo;
    for i = 2:n
        x(i,:) = x(i-1,:) + (hi-x(i-1,:))./((n-i+1)^phi);
    end

end
function par = update_struct(par,names,vals)
    
    if iscell(vals)
        for i = 1:numel(names)
            par.(names{i}) = vals{i};
        end    
    else
        for i = 1:numel(names)
            par.(names{i}) = vals(i);
        end
    end
end
function covAB = cov(A,B)
   
    covmat = nancov(A,B);
    covAB = covmat(1,2);
    
end

end
end

