function [minimums, maximums] = FVA_own(model, rxns, factor)

if nargin<2
    rxns = model.rxns;
end

if nargin<3
    factor = 1;
end

fba = optimizeCbModel(model);
minimums = zeros(size(rxns));
maximums = zeros(size(rxns));
for i = 1:length(rxns)
    modelAux = changeObjective(model,rxns(i));
    modelAux = changeRxnBounds(modelAux,model.rxns(model.c==1) , fba.f*factor,'l');
    mi = optimizeCbModel(modelAux,'min');
    ma = optimizeCbModel(modelAux,'max');
    minimums(i) = mi.f;
    maximums(i) = ma.f;
end


end