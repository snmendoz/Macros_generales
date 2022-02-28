function [posExc, excRxn] = getExchangeRxnFromMetabolite(model, metabolite)
excRxn = {};
pos_ExcRxns = find(findExcRxns(model));
[~, posRxnsByMet, ~, ~, ~] = getRxnsFromMets(model, metabolite);
posExc = intersect(pos_ExcRxns,posRxnsByMet);
if ~isempty(posExc)
    excRxn = model.rxns{posExc};
end
end