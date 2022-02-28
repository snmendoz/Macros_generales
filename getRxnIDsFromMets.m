function [reactionsByMet,posReactionsByMet, allReactions, posAllReactions] = getRxnIDsFromMets(model, mets)

if ischar(mets)
    pos = find(strcmp(model.mets, mets));
elseif iscell(mets)
    pos = cell2mat(arrayfun(@(x)find(strcmp(x,model.mets)),mets,'UniformOutput',false))';
else
    pos = mets;
end

n_mets = length(pos);

if iscell(mets)
    reactionsByMet = cell(size(mets));
    posReactionsByMet = cell(size(mets));
    allReactions = {};
    posAllReactions = {};
    for i = 1:n_mets
        posRxns = find(model.S(pos(i), :));
        posReactionsByMet{i} = posRxns;
        reactionsByMet{i} = model.rxns(posRxns);
        posAllReactions = union(posAllReactions, posReactionsByMet{i});
        allReactions = union(allReactions,reactionsByMet{i});
    end
else
    posRxns = find(model.S(pos, :));
    posReactionsByMet = posRxns;
    posAllReactions = posReactionsByMet;
    reactionsByMet = model.rxns(posRxns);
    allReactions = reactionsByMet;
end

end