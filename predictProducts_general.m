function productExchangeRxns = predictProducts_general(model, exchangeRxnDetectionMethod, option, growthRateConstraintFraction, tol)

if nargin < 5
    tol = 0;
end
% option 1 returns only products which MUST be produced
% option 2 returns all the metabolites that can eventually be produced
if exchangeRxnDetectionMethod ==1
    exch = model.rxns(findExcRxns(model));
elseif exchangeRxnDetectionMethod ==2
    exch = model.rxns(findExcRxnsWithIDs(model));
end

fba_basal = optimizeCbModel(model);

minimums = zeros(size(exch));
maximums = zeros(size(exch));
    
for i = 1:length(exch)
    modelAux = changeObjective(model,exch(i));
    modelAux = changeRxnBounds(modelAux,model.rxns(find(model.c)) , fba_basal.f*growthRateConstraintFraction,'l');
    mi = optimizeCbModel(modelAux,'min');
    ma = optimizeCbModel(modelAux,'max');
    minimums(i) = mi.f;
    maximums(i) = ma.f;
end

switch option    
    case 1
        productExchangeRxns = exch(intersect(find(maximums(:,1)>tol), find(minimums(:,1)>tol)));
    case 2
        productExchangeRxns = exch(find(maximums(:,1)>tol));
end

end