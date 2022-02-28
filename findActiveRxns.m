function [rxns, posRxns] = findActiveRxns(model, fluxDistribution)

posRxns = 1:length(model.rxns);

rxns = {};
for i = 1:length(posRxns)
    if model.lb(posRxns(i))~=0
        if model.lb(posRxns(i)) == fluxDistribution(posRxns(i))
            rxns = union(rxns,model.rxns(posRxns(i)));
        end
    end
end

posRxns = getPosOfElementsInArray(rxns, model.rxns);

end