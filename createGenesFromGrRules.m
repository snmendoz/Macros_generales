function model = createGenesFromGrRules(model)
grRules = model.grRules;
nonEmptyPos = find(cellfun(@isempty, grRules)==0);
grRules = grRules(nonEmptyPos);
fields = fieldnames(model);
geneFields = fields(find(~cellfun(@isempty, regexp(fields, '^gene'))));
geneFields = geneFields(cellfun(@(x) iscell(model.(x)), geneFields));
geneFields = setdiff(geneFields, 'genes');
if isfield(model, 'proteins')
    geneFields = union(geneFields, 'proteins');
end

involvedGenes = cellfun(@unique,cellfun(@splitString,regexprep(grRules,'(?<=\ )or(?=\ )|(?<=\ )and(?=\ )|\(|\)',' '),'UniformOutput',0),'UniformOutput',0);
genes = {};
for i=1:length(involvedGenes) 
    genes = union(genes,involvedGenes{i});
end

if size(genes,2)>size(genes,1)
    genes = genes';
end

if ~isfield(model, 'genes')
    model.genes = {};
end

inter = intersect(genes, model.genes);
for i = 1:length(geneFields)
    newField = repmat({''}, length(genes),1);
    newField(getPosOfElementsInArray(inter,genes)) = model.(geneFields{i})(getPosOfElementsInArray(inter,model.genes));
    model.(geneFields{i}) = newField;
end

model.genes = genes;
model = createRulesFromgrRules(model);
model = create_rxnGeneMat(model);

end