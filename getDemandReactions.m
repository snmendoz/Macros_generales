function DM_rxns = getDemandReactions(model)

DM_rxns = model.rxns(startsWith(model.rxns,'DM_'));

end