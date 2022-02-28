function [model, rxnDiff ,fba] = forceBiomassProduction(draft, original)

rxnDiff = setdiff(original.rxns, draft.rxns);
model = draft;
for i = 1:length(rxnDiff)
    eq = getRxn_cobraFormat(original, find(strcmp(original.rxns, rxnDiff{i})));
    model = addReaction(model, rxnDiff{i}, 'reactionFormula', eq{1});
end

% [int, pos1, pos2] = intersect(original.rxns, model.rxns);
% model.lb(pos2) = original.lb(pos1);
% model.ub(pos2) = original.ub(pos1);
% for i = 1:length(int)
% 
%     [sameRxnID, sameRxnFormula, writtenInSameSense, sameReversibility] = compareRxns(model, pos2(i), original, pos1(i));
%     if ~sameRxnFormula || ~writtenInSameSense
%         disp('')
%     end
% end
    

% [rxnDiff2, pos] = setdiff(draft.rxns, original.rxns);
% a = find(draft.ub(pos)<0);
% b = find(draft.lb(pos)>0);

fba = optimizeCbModel(model);

end