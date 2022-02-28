function relateGenesToRxn3(models, idsType, geneMatrix, rxnMatrix, rxnIDs, rxnMatrix2, rxnIDs2, metOtherIDs, metMNXIDs, rxnOtherIDs, rxnMNXIDs, species)

manualModel = models{1};
n_genes = length(manualModel.genes); 
n_rxns = length(manualModel.rxns);

rxnMatrix3 = rxnMatrix;
rxnIDs3 = rxnIDs;
for i = 1:n_rxns
    for j = 3:length(models)-1
        if rxnMatrix2(i,j)==1
            rxnMatrix3(i,j)=1;
            rxnIDs3(i,j)=rxnIDs2(i,j);
        end
    end
end
  
gene_rxns_Matrix = cell(size(geneMatrix));
gene_rxns_Matrix_specific = cell(size(geneMatrix));

percentagesNotRecovered = cell(size(geneMatrix));
percentagesAdditional = cell(size(geneMatrix));

recoveredRxns_perGene_perModel = cell(length(models)-1,1);
additionalRxns_perGene_perModel = cell(length(models)-1,1);

for i = 2:length(models) 
    pos_genes_i = find(geneMatrix(1:n_genes,i));
    genes_i = manualModel.genes(find(geneMatrix(1:n_genes,i)));
    recoveredRxns_perGene = 100*ones(size(genes_i));
    additionalRxns_perGene = zeros(size(genes_i));
    
    model_i = models{i};
    if isfield(model_i, 'grRules') && ~isfield(model_i, 'rules')
        model_i = createRulesFromgrRules(model_i);
    end
    model_i = creategrRulesField(model_i);

    if i ==1
       disp('') 
    end
    
    for j = 1:n_genes
        fprintf(['i:' num2str(i) ' j:' num2str(j) ' \n'])
        relatedRxns_manual = manualModel.rxns(find(~cellfun(@isempty, strfind(manualModel.grRules, manualModel.genes{j}))));
        relatedRxns_i = model_i.rxns(find(~cellfun(@isempty, strfind(model_i.grRules, manualModel.genes{j}))));

        if isempty(relatedRxns_manual)
           disp('') 
        end
        
        logical_rxnWasRecovered = zeros(size(relatedRxns_manual));
        logical_rxnIsAdditional = ones(size(relatedRxns_i));
        pos_match = zeros(size(relatedRxns_manual));
        n_recovered = 0;
        for k = 1:length(relatedRxns_manual)
            pos = find(strcmp(rxnIDs3(1:n_rxns,2) ,relatedRxns_manual{k}));
            rxn = rxnIDs3(pos, i+1);
            if ~isempty(rxn) 
                [is, posMatch] = ismember(rxn,relatedRxns_i);
                if is
                    n_recovered = n_recovered+1;
                    logical_rxnWasRecovered(k) = 1;
                    pos_match(k) = posMatch;
                    logical_rxnIsAdditional(posMatch) = 0;
                end
            end
        end
        n_additional = length(find(logical_rxnIsAdditional));
        
        st1 = '';
        st2 = '';
        if ~isempty(find(logical_rxnWasRecovered==0, 1))
            notRecovered = relatedRxns_manual(logical_rxnWasRecovered==0);
            eqs_notRecovered = getRxn_cobraFormat(manualModel, notRecovered);
            st_1 = strcat(notRecovered, ':' ,eqs_notRecovered);
            st1 = strjoin (st_1,', ');
            st1 = strcat('notRecovered: ', st1);
            percentagesNotRecovered{j,i} = 100*n_recovered/length(relatedRxns_manual);
            if ismember(j, pos_genes_i)
                recoveredRxns_perGene(pos_genes_i==j) = percentagesNotRecovered{j,i};
            end
        end
        if n_additional > 0
            additionalRxns = relatedRxns_i(logical_rxnIsAdditional==1);
            eqs_additional = getRxn_cobraFormat(model_i, additionalRxns);
            st_1 = strcat(additionalRxns, ':' ,eqs_additional);
            st2 = strjoin (st_1,', ');
            st2 = strcat('. Additional: ', st2);
            percentagesAdditional{j,i} = 100*n_additional/length(relatedRxns_manual);
            if ismember(j, pos_genes_i)
                additionalRxns_perGene(pos_genes_i==j) = percentagesAdditional{j,i};
            end
        end
        gene_rxns_Matrix{j,i} = strcat(st1, st2);
        if ismember(j, pos_genes_i)
            gene_rxns_Matrix_specific{j,i} = strcat(st1, st2);
        end
    end
    recoveredRxns_perGene_perModel{i-1} = recoveredRxns_perGene;
    additionalRxns_perGene_perModel{i-1} = additionalRxns_perGene;
end

averageRecoveredRxns_perGene_perModel = cellfun(@(x) mean(x),recoveredRxns_perGene_perModel);
averageAdditionalRxns_perGene_perModel = cellfun(@(x) mean(x),additionalRxns_perGene_perModel);

save(['gene_rxns_Matrix_' species], 'gene_rxns_Matrix');
save(['percentagesNotRecovered_' species], 'percentagesNotRecovered');
save(['percentagesAdditional_' species], 'percentagesAdditional');

if exist(['geneRxnRelations_' species '.xlsx'],'file')==2
    delete(['geneRxnRelations_' species '.xlsx']);
end

xlswrite(['geneRxnRelations_' species '.xlsx'],averageRecoveredRxns_perGene_perModel,'avRec');
xlswrite(['geneRxnRelations_' species '.xlsx'],averageAdditionalRxns_perGene_perModel,'avAdd');
xlswrite(['geneRxnRelations_' species '.xlsx'],gene_rxns_Matrix,'gene_rxns_Matrix');
xlswrite(['geneRxnRelations_' species '.xlsx'],gene_rxns_Matrix_specific,'gene_rxns_Matrix_sp');
xlswrite(['geneRxnRelations_' species '.xlsx'],percentagesNotRecovered,'percentagesNotRecovered');
xlswrite(['geneRxnRelations_' species '.xlsx'],percentagesAdditional,'percentagesAdditional');


end