function model = changeCoefficient(model, rxn, met, coef)

model.S(strcmp(model.mets,met)==1, strcmp(model.rxns,rxn)==1) = coef;

end