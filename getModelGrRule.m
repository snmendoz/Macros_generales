function grRule = getModelGrRule(model,rxn)

grRule = model.grRules{getPosOfElementsInArray({rxn}, model.rxns)};

end