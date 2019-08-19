function [] = compile_mex(name,threads,onedim,intel)
           
    % a. clear existing mex-file in userpath
    dirpath = sprintf('%s/cfuncs',userpath);
    if isdir(dirpath)
        rmdir(dirpath,'s')
    end
    copyfile(sprintf('%s/cfuncs',pwd),sprintf('%s/cfuncs',userpath))
    if onedim == 1
        delete(sprintf('%s\\%s_onedim.mexw64',userpath,name))        
    end
    
    % b. build string
    str = sprintf('mex -largeArrayDims -outdir %s %s/cfuncs/%s.cpp -DMAXTHREADS=%d',userpath,userpath,name,threads);
    if strcmp(name,'mex_solve') && onedim == 1
        str = sprintf('mex -largeArrayDims -outdir %s -output mex_solve_onedim %s/cfuncs/%s.cpp -DONEDIM=1 -DMAXTHREADS=%d',...
            userpath,userpath,name,threads);    
    end
    
        % flags
        if intel == 0
            str = sprintf('%s %s',str,' CXXFLAGS="$CXXFLAGS -std=c++11 -O3 -Wall -fopenmp -ffast-math"');
        else
            str = sprintf('%s %s',str,' COMPFLAGS=''$COMPFLAGS /o3 /openmp /arch:CORE-AVX512''');          
        end
        
        % libgomp (for OpenMP)
        if intel == 0
            str = sprintf('%s %s',str,' C:/ProgramData/MATLAB/SupportPackages/R2018b/3P.instrset/mingw_w64.instrset/lib/gcc/x86_64-w64-mingw32/6.3.0/libgomp.a'); 
        end
        
    % c. evaluate
    eval(str);
    fprintf('\n');
    
end