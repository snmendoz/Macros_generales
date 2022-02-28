function [bigg, removed] = removeDuplicatedReactionEquationsInBiGG(bigg)
if nargin<1
    load('D:\Dropbox\Databases\BIGG\bigg_85.mat')
end

[~, redundantRxnPos, ~] = identifyDuplicatedReactionsInBiGG(bigg);

removed = {};
for i = 1:length(redundantRxnPos)
    rxns_i = bigg.rxns(redundantRxnPos{i});
    contains_copy = cellfun(@(x) contains(x, 'copy'), rxns_i);
    if any(contains_copy)
        if length(find(contains_copy==0))==1
            removed = union(removed, rxns_i(contains_copy));
        elseif length(find(contains_copy==0))==0     
            rxns_i_sorted = sort(rxns_i);
            removed = union(removed, rxns_i_sorted(2:end));
        else
            pos1 = find(contains_copy); pos1 = pos1(1);
            copy1 = rxns_i(pos1);
            backbone = regexprep(copy1,'_copy.*','');
            diff = setdiff(rxns_i,backbone);
            removed = union(removed, diff);
        end
    else
        removed = union(removed, rxns_i(2:end));
    end
end

bigg = removeRxns(bigg, removed);
save('D:\Dropbox\Databases\BIGG\bigg_85_uniqueRxnEqs.mat', 'bigg')
end