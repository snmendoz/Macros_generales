function [metIDs, allMetIDs, allMetNames, allMetPos] = getMetaboliteIDsFromRxns(model, reactions)

pos_reactions = cell2mat(arrayfun(@(x)find(strcmp(x,model.rxns)),reactions,'UniformOutput',false))';

allMetIDs = {};
metIDs = cell(size(reactions));

for i = 1:length(reactions)
    mets = model.mets(model.S(:,pos_reactions(i))~=0);
    metIDs{i} = mets;
    allMetIDs = union(allMetIDs, mets);
end

allMetPos = cell2mat(arrayfun(@(x)find(strcmp(x,model.mets)),allMetIDs,'UniformOutput',false))';
allMetNames = model.metNames(allMetPos);
end