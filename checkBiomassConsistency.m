function total = checkBiomassConsistency(model)

biomass = model.rxns{find(model.c)};

mets = model.mets(find(model.S(:, getPosOfElementsInArray({biomass}, model.rxns))));
mws = cellfun(@(x) calculateMolecularWeight(model, x), mets);
coefs = full(model.S(find(model.S(:, getPosOfElementsInArray({biomass}, model.rxns))), getPosOfElementsInArray({biomass}, model.rxns)));

total = -sum(coefs.*mws)/1000;

end