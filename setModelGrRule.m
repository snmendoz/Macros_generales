function model = setModelGrRule(model, rxn, newRule)

model.grRules{getPosOfElementsInArray({rxn}, model.rxns)} = newRule;

end