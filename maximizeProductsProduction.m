function [products,isProducedByNetwork, productsNotProduced] = maximizeProductsProduction(model, rxn)

products = getProducts(model, rxn);
isProducedByNetwork = zeros(size(products));
howMuch = zeros(size(products));
for i = 1:length(products)
    modelAux = model;
    modelAux = addReaction(modelAux,['ownDM_' products{i}],...
        'reactionFormula',[products{i}  ' -> ']);
    modelAux = changeObjective(modelAux, ['ownDM_' products{i}]);
    fba_i = optimizeCbModel(modelAux);
    if fba_i.f>0
        isProducedByNetwork(i) = 1;
        howMuch(i) = fba_i.f;
    end
end
productsNotProduced = products(find(isProducedByNetwork==0));

end