function limitingNutrients = findLimitingNutrient(model,type)

if nargin < 3
    type ='relative';
end

perturbance_parameter_absolute = 0.01;
perturbance_parameter_relative = 1.01;

fba = optimizeCbModel(model);
[rxns, posRxns] = findActiveExchangeRxns(model, fba.x);
limitingNutrients = {};

for i = 1:length(rxns)
    modelAux = model;
    if strcmp(type, 'relative')
        modelAux = changeRxnBounds(modelAux,rxns{i},modelAux.lb(posRxns(i))*perturbance_parameter_relative,'l');
    elseif strcmp(type, 'absolute')
        modelAux = changeRxnBounds(modelAux,rxns{i},modelAux.lb(posRxns(i))*-perturbance_parameter_absolute,'l');
    end
    fbaAux = optimizeCbModel(modelAux);
    if fbaAux.f>fba.f+10^-8 || fbaAux.f>fba.f*1.001
        limitingNutrients = union(limitingNutrients, rxns(i));
    end
end

end