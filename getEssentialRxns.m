function essentials = getEssentialRxns(model, threshold,rxnList)

fba = optimizeCbModel(model);
if nargin < 2
   threshold = fba.f*0.1; 
end

if (nargin < 3)
    rxnList = model.rxns;
else
    if (isempty(rxnList))
        rxnList = model.rxns;
    end
end

[grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(model, 'FBA', rxnList, 'false');


posEx = find(cellfun(@isempty, strfind(model.rxns, 'EX_'))==0);
n_essentials = 0;
essentials = {};

for i = 1:length(posEx)

   model2 = changeRxnBounds(model,model.rxns(posEx(i)),0,'l') ;
   fba2 = optimizeCbModel(model2);
   if fba2.f < threshold
       n_essentials = n_essentials + 1;
       essentials = union(essentials, model.rxns(posEx(i)));
   end
end

end