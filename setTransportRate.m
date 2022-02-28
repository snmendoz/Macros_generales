function model = setTransportRate(model, rxns, type, value, bound)

for i = 1:length(rxns)
    disp(i)
    [met, extracellular_in_left] = getTransportedMetabolite(model,rxns{i});
    if ~isempty(met)
        if strcmp(type, 'uptake')
            if extracellular_in_left
                if strcmp(bound,'max')
                    model = changeRxnBounds(model, rxns{i}, value, 'u');
                elseif strcmp(bound,'min')
                    model = changeRxnBounds(model, rxns{i}, value, 'l');
                elseif strcmp(bound,'equal')
                    model = changeRxnBounds(model, rxns{i}, value, 'b');
                end
            else
                if strcmp(bound,'max')
                    model = changeRxnBounds(model, rxns{i}, -value, 'l');
                elseif strcmp(bound,'min')
                    model = changeRxnBounds(model, rxns{i}, -value, 'u');
                elseif strcmp(bound,'equal')
                    model = changeRxnBounds(model, rxns{i}, -value, 'b');
                end
            end
        elseif strcmp(type, 'secretion')
            if extracellular_in_left
                if strcmp(bound,'max')
                    model = changeRxnBounds(model, rxns{i}, -value, 'l');
                elseif strcmp(bound,'min')
                    model = changeRxnBounds(model, rxns{i}, -value, 'u');
                elseif strcmp(bound,'equal')
                    model = changeRxnBounds(model, rxns{i}, -value, 'b');
                end
            else
                if strcmp(bound,'max')
                    model = changeRxnBounds(model, rxns{i}, value, 'u');
                elseif strcmp(bound,'min')
                    model = changeRxnBounds(model, rxns{i}, value, 'l');
                elseif strcmp(bound,'equal')
                    model = changeRxnBounds(model, rxns{i}, value, 'b');
                end
            end
        end
        
    end
end

end

