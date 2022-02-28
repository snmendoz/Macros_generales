function [x,y,x_adicionales,y_adicionales] = runBensolve(runID, model, obj, outputFileName, del)

current = pwd;
bensolveFullPath = which('bensolve.exe'); 
bensolveDirectory = fileparts(bensolveFullPath);

if ~isempty(bensolveDirectory); 
    cd(bensolveDirectory)
    if ~isdir(runID); mkdir(runID); end; cd(runID)
    runPath = pwd;
    %export inputs to bensolve
    if ~isdir('Inputs'); mkdir('Inputs'); end; cd('Inputs')
    if del
        delete([outputFileName '_img*'])
        delete([outputFileName '_adj*'])
        delete([outputFileName '_inc*'])
        delete([outputFileName '_pre*'])
        delete([outputFileName '.vlp'])
        delete([outputFileName '.log'])
        delete([outputFileName '._c*'])
        delete([outputFileName '_c.sol'])
    end
    cobra2vlp(model, obj, outputFileName)
    if ~exist([outputFileName '_img_p.sol'])
        cd(bensolveDirectory)
        run = system(['bensolve.exe ' runID '/Inputs/' outputFileName '.vlp -s']);
        if run ~= 0 && run ~= 1
            disp(['run this: "./bensolve.exe ' runID '/Inputs/' outputFileName '.vlp -s" and then press enter']  )
            pause;
        elseif run ==0
            disp('Bensolve was executed correctly')
            dir = [bensolveDirectory '\' runID '\Inputs'];
            [x,y,x_adicionales,y_adicionales] = ExtraerPuntosPareto(dir,[outputFileName '_img_p.sol']);
        else 
            disp('The problem is infeasible')
            x = []; y = []; x_adicionales = []; y_adicionales = [];
        end
    else
        dir = [bensolveDirectory '\' runID '\Inputs'];
        [x, y, x_adicionales,y_adicionales] = ExtraerPuntosPareto(dir,[outputFileName '_img_p.sol']);
    end    
else
end

cd(current);




end