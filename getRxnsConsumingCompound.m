function [rxnsConsumingATP, fluxRxnsConsumingATP] = getRxnsConsumingCompound(model,compound, FBAsolution)

if allMetsInCBMPYFormat(model.mets)
    atp_molecule = transformMetsToCBMPYFormat(compound);
else
    atp_molecule = transformMetsToCOBRAFormat(compound);

end
atp_molecule =atp_molecule{1};

posATP = getPosOfElementsInArray({atp_molecule}, model.mets);
rxns = getRxnsFromMets(model, atp_molecule);
posRxns = getPosOfElementsInArray(rxns, model.rxns);
isReactant = arrayfun(@(x) full(model.S(posATP,posRxns(x)))<0, 1:length(posRxns));
isFluxPositiveDirection = FBAsolution(posRxns)>0;
isFluxNegativeDirection = FBAsolution(posRxns)<0;

posRxnsProducingATP = union(intersect(find(isReactant), find(isFluxPositiveDirection)), intersect(find(~isReactant), find(isFluxNegativeDirection)));
rxnsConsumingATP = rxns(posRxnsProducingATP);
fluxRxnsConsumingATP = FBAsolution(posRxns(posRxnsProducingATP));

end