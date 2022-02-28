function model = setTransportUptake(model, rxns, value, bound)

for i = 1:length(rxns)
    disp(i)
    [met, extracellular_in_left] = getTransportedMetabolite(model,rxns{i});
    if ~isempty(met)
        if extracellular_in_left
            if strcmp(bound,'max')
                model = changeRxnBounds(model, rxns{i}, value, 'u');
            elseif strcmp(bound,'min')
                model = changeRxnBounds(model, rxns{i}, value, 'l');
            end
            
        else
            
            if strcmp(bound,'max')
                model = changeRxnBounds(model, rxns{i}, -value, 'l');
            elseif strcmp(bound,'min')
                model = changeRxnBounds(model, rxns{i}, -value, 'u');
            end
            
        end
    end
end

end

