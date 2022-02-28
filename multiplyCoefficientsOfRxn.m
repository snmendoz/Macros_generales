function model = multiplyCoefficientsOfRxn(model, rxn, factor)

mets = getMetaboliteIDsFromRxns(model, {rxn});
mets = mets{1};
for i = 1:length(mets)
    newCoeff = full(model.S(getPosOfElementsInArray(mets(i),model.mets),getPosOfElementsInArray({rxn},model.rxns)))*factor;
    model = changeCoefficient(model, rxn,mets{i},newCoeff);
end


end