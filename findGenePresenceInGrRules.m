function pos_present = findGenePresenceInGrRules(model, gene)

grRules = model.grRules;
involvedGenes = cellfun(@unique,cellfun(@splitString,regexprep(grRules,{'\(|\)','\ or\ |\ and\ |\ &\ |\ \|\ '},{'',' '}),'UniformOutput',0),'UniformOutput',0);
pos_present = find(cellfun(@(x) ismember(gene,x), involvedGenes));


end