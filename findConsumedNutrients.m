function consumed = findConsumedNutrients(model)

fba = optimizeCbModel(model);
pos_exchange = find(findExcRxns(model));
consumed = model.rxns(intersect(pos_exchange, find(fba.x<0)));

end