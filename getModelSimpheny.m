function model = getModelSimpheny(fileName)

[rxnIDs, rxnNames, equations, geneAssociations] = getRxnsFromSimpheny(fileName);

load('D:\Dropbox\Databases\BIGG\e_coli_core');
model = e_coli_core;
model.rxns = strcat('OWN_', model.rxns);
rxns = model.rxns;

for i = 1:length(rxnIDs)
    model = addReaction(model, rxnIDs{i}, 'reactionFormula', equations{i}, 'geneRule', geneAssociations{i}); 
end
model = removeRxns(model, rxns);
model = rmfield(model, 'rev');
model = createGenesFromGrRules(model);
model = createRulesFromgrRules(model);

end