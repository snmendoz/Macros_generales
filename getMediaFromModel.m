function [media, summary] = getMediaFromModel(model,useIds)

if useIds
    [posExcRxn, excRxns] = findExcRxnsWithIDs(model);
else
    posExcRxn = find(findExcRxns(model));
    excRxns = modle.rxns(posExcRxn);
end

lb_exc = model.lb(posExcRxn);
ub_exc = model.ub(posExcRxn);

pos_nutrients = find(lb_exc<0);

media.reactions = excRxns(pos_nutrients);
media.lb = lb_exc(pos_nutrients);
media.ub = ub_exc(pos_nutrients);

summary = [media.reactions num2cell(media.lb) num2cell(media.ub)];

end