function mets = getMetsOfDemandReactions(model)

mets = {};
DM_rxns = model.rxns(startsWith(model.rxns,'DM_'));
if ~isempty(DM_rxns)
    [~, mets,~, ~] = getMetaboliteIDsFromRxns(model, DM_rxns);
end

end