function [essentials,grRateKO,pos_essentials] = findEssentialRxns(model, threshold,rxnList)

fba = optimizeCbModel(model);
if nargin < 2 || isempty(threshold)
   threshold = fba.f*0.1; 
end

if (nargin < 3) || isempty(threshold)
    rxnList = model.rxns;
end

[grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(model, 'FBA', rxnList, 'false');

essentials = grRateKO<threshold;
pos_essentials = getPosOfElementsInArray(rxnList(essentials), model.rxns); 

end