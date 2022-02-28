function [rxnsConsumingCompound, fluxRxnsConsumingCompound] = getRxnsConsumingCompound_FVA(model,compound, FBAsolution)

if allMetsInCBMPYFormat(model.mets)
    compound = transformMetsToCBMPYFormat(compound);
else
    compound = transformMetsToCOBRAFormat(compound);

end
compound =compound{1};

posCompound = getPosOfElementsInArray({compound}, model.mets);
rxns = getRxnsFromMets(model, compound);

[minimums, maximums] = FVA_own(model, rxns);

posRxns = getPosOfElementsInArray(rxns, model.rxns);
isReactant = arrayfun(@(x) full(model.S(posCompound,posRxns(x)))<0, 1:length(posRxns));
canReactionBeForward = maximums>1e-9;
canReactionBeBackward = minimums<-1e-9;

posRxnsConsumingCompound = union(intersect(find(isReactant), find(canReactionBeForward)), intersect(find(~isReactant), find(canReactionBeBackward)));
rxnsConsumingCompound = rxns(posRxnsConsumingCompound);
fluxRxnsConsumingCompound = FBAsolution(posRxns(posRxnsConsumingCompound));

end