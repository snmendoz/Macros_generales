function [rxns, pos] = getInterExchangeReactions_steadyCom(communityModel)
%just those which have [e] and [u]
has = cellfun(@(x) rxnHasCompartment(communityModel, x, {'u','e'}), communityModel.rxns);
pos = find(has);
rxns = communityModel.rxns(pos);

end