function connectivity = getMetaboliteConnectivityForMet(model, met)

connectivity = length(find(model.S(getPosOfElementsInArray({met},model.mets),:)));
end