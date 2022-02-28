function total = getCHONSfromUptake(model,excRxn,flux)

total = 0;
for i = 1:length(excRxn)
    excRxn_i = excRxn(i);
    pos_met = getMetFromExcRxn(model, getPosOfElementsInArray(excRxn_i, model.rxns));
    metFormula = model.metFormulas{pos_met};
    [Ematrix, elements] = getElementalComposition(metFormula, {'C','H','O','N','S'});
    total = total+Ematrix*flux(i);
end

end