function [reactants,isProducedByNetwork, reactantsNotProduced] = maximizeReactantsProduction(model, rxn)

reactants = getReactants(model, rxn);
isProducedByNetwork = zeros(size(reactants));
howMuch = zeros(size(reactants));
for i = 1:length(reactants)
    modelAux = model;
    modelAux = addReaction(modelAux,['ownDM_' reactants{i}],...
        'reactionFormula',[reactants{i}  ' -> ']);
    modelAux = changeObjective(modelAux, ['ownDM_' reactants{i}]);
    fba_i = optimizeCbModel(modelAux);
    if fba_i.f>0
        isProducedByNetwork(i) = 1;
        howMuch(i) = fba_i.f;
    end
end
reactantsNotProduced = reactants(find(isProducedByNetwork==0));

end