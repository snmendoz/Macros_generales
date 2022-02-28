function [rxnsProducingATP, posRxnsProducingATP, fluxRxnsProducingATP] = getRxnsProducingATP(model, FBAsolution)

if allMetsInCBMPYFormat(model.mets)
    atp_molecule = 'atp_c';
else
    atp_molecule = 'atp[c]';

end

posATP = getPosOfElementsInArray({atp_molecule}, model.mets);
rxns = getRxnsFromMets(model, atp_molecule);
posRxns = getPosOfElementsInArray(rxns, model.rxns);
isReactant = arrayfun(@(x) full(model.S(posATP,posRxns(x)))<0, 1:length(posRxns));
isFluxPositiveDirection = FBAsolution(posRxns)>0;
isFluxNegativeDirection = FBAsolution(posRxns)<0;

posRxnsProducingATP = union(intersect(find(isReactant), find(isFluxNegativeDirection)), intersect(find(~isReactant), find(isFluxPositiveDirection)));
rxnsProducingATP = rxns(posRxnsProducingATP);
fluxRxnsProducingATP = FBAsolution(posRxns(posRxnsProducingATP));

end