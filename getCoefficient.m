function coefficient = getCoefficient(model, rxn, met)

coefficient = full(model.S(getPosOfElementsInArray({met}, model.mets), getPosOfElementsInArray({rxn}, model.rxns)));

end