function [rxns, posRxns] = findActiveExchangeRxns(model, fluxDistribution)

posExchRxns = find(findExcRxns(model));

rxns = {};
for i = 1:length(posExchRxns)
    if model.lb(posExchRxns(i))~=0
        if model.lb(posExchRxns(i)) == fluxDistribution(posExchRxns(i))
            rxns = union(rxns,model.rxns(posExchRxns(i)));
        end
    end
end

posRxns = getPosOfElementsInArray(rxns, model.rxns);

end