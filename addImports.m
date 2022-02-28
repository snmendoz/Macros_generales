function model = addImports(model, seeds)

for i = 1:length(seeds)
    rxnID  = ['import_' seeds{i}];
    formula = [seeds{i} '_e' ' <==> ' seeds{i} '_c'];
    model = addReaction(model,rxnID, 'reactionFormula', formula);
    
    rxnID = ['EX_' seeds{i} '_e'];
    model = addReaction(model, rxnID, 'metaboliteList', {[seeds{i} '_e']}, 'stoichCoeffList', -1);
    
end

end