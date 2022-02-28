function saveFluxValues(variableName, values)

if exist([variableName '.mat'],'file')==2
    load([variableName '.mat']);
    fluxes = [fluxes, values];
    save([variableName '.mat'], 'fluxes');

else
    fluxes = values;
    save([variableName '.mat'], 'fluxes');
end


end