function coef = getStoichiometricCoefficient(model, rxn, met)

coef = full(model.S(getPosOfElementsInArray({met}, model.mets), getPosOfElementsInArray({rxn}, model.rxns)));

end