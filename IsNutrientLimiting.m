function is = IsNutrientLimiting(model, nutrient,type)

if nargin < 3
    type ='relative';
end

perturbance_parameter_absolute = 0.01;
perturbance_parameter_relative = 1.01;

is = 0;
fba = optimizeCbModel(model);
modelAux = model;
posRxn = getPosOfElementsInArray({nutrient}, modelAux.rxns);
if strcmp(type, 'relative')
    modelAux = changeRxnBounds(modelAux,nutrient,modelAux.lb(posRxn)*perturbance_parameter_relative,'l');
elseif strcmp(type, 'absolute')
    modelAux = changeRxnBounds(modelAux,nutrient,modelAux.lb(posRxn)-perturbance_parameter_absolute,'l');
end
fbaAux = optimizeCbModel(modelAux);
if fbaAux.f>fba.f+10^-8 || fbaAux.f>fba.f*1.001
    is = 1;
end

end