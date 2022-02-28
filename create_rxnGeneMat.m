function model = create_rxnGeneMat(model)


rxnGeneMat=zeros(length(model.rxns),length(model.genes));
genes_Involucrados_por_Rxn=regexp(regexprep(model.grRules,'\or|and|\(|\)',''),'\  ','split');

for i=1:length(model.rxns) 
    
    genesInvolucrados=genes_Involucrados_por_Rxn{i};
    [~,indices]=intersect(model.genes,genesInvolucrados);
    for j=1:length(indices)
        rxnGeneMat(i,indices(j))=1;
    end   
end

model = setfield(model,'rxnGeneMat',rxnGeneMat);
