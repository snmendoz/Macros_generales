function bounds = getBoundsForEquivalentDemandReactins(model,rxnList, mets, fba)

bounds = zeros(size(rxnList));
for i = 1:length(rxnList)
    disp(i)
    internal_rxn_i = rxnList{i};

    coef = getStoichiometricCoefficient(model, internal_rxn_i, mets{i});
    flux_aa = fba.x(getPosOfElementsInArray({internal_rxn_i}, model.rxns));
    
    bounds(i) = coef*flux_aa;
end