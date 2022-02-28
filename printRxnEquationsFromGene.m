function printRxnEquationsFromGene(model, gene, oldLocusTag)

if nargin <3; oldLocusTag = 0; end;

if oldLocusTag
   posOld = getPosOfElementsInArray({gene}, model.geneOldLocusTag);
   if ~isempty(posOld)
       gene = model.genes{posOld};
   else
       return;
   end
end

[~, pos] = search(gene, model.grRules);

allEquations = getRxn_cobraFormat(model, pos);
allEquations = strcat(model.rxns(pos),' : ',allEquations);
for i = 1:size(allEquations,1)
    fprintf([allEquations{i} '\n']);
end
end