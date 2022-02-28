function relateGenesToRxns(models, idsType, geneMatrix, rxnMatrix, rxnIDs, rxnMatrix2, rxnIDs2, metOtherIDs, metMNXIDs, rxnOtherIDs, rxnMNXIDs, species)

manualModel = models{1};
n_genes = length(manualModel.genes); 
n_rxns = length(manualModel.rxns);
  
gene_rxns_Matrix = cell(size(geneMatrix));
percentagesNotRecovered = cell(size(geneMatrix));
percentagesAdditional = cell(size(geneMatrix));


for i = 2:length(models) 
    genes_i = manualModel.genes(find(geneMatrix(1:n_genes,i)));
    recoveredRxns_perGene = 100*ones(size(genes_i));
    additionalRxns_perGene = zeros(size(genes_i));

    model_i = models{i};
    if isfield(model_i, 'grRules') && ~isfield(model_i, 'rules')
        model_i = createRulesFromgrRules(model_i);
    end
    model_i = creategrRulesField(model_i);
    if strcmp(idsType{i},'bigg')
        
    else
        model_i = removeDuplicatedMetabolitesFromModel(model_i);
        model_i = removeEmptyRxns(model_i);
        model_i = translateModelToTargetLanguage(model_i, idsType{i}, 'bigg', manualModel, metOtherIDs, metMNXIDs, rxnOtherIDs, rxnMNXIDs, 0);
    end

    for j = 1:n_genes
        fprintf(['i:' num2str(i) ' j:' num2str(j) ' \n'])
        relatedRxns_manual = manualModel.rxns(find(~cellfun(@isempty, strfind(manualModel.grRules, manualModel.genes{j}))));
        relatedRxns_i = model_i.rxns(find(~cellfun(@isempty, strfind(model_i.grRules, manualModel.genes{j}))));

        recoveredRxns_perGene(j) = 100*length(find(ismember(relatedRxns_manual, relatedRxns_i)))/length(relatedRxns_manual);
        additionalRxns_perGene(j) = 100*length(setdiff(relatedRxns_i, relatedRxns_manual))/length(relatedRxns_manual);
        
        st1 = '';
        st2 = '';
        if recoveredRxns_perGene(j) < 100
            notRecovered = setdiff(relatedRxns_manual, relatedRxns_i);
            eqs_notRecovered = getRxn_cobraFormat(manualModel, notRecovered);
            st_1 = strcat(notRecovered, ':' ,eqs_notRecovered);
            st1 = strjoin (st_1,', ');
            st1 = strcat('notRecovered: ', st1);
            percentagesNotRecovered{j,i} = recoveredRxns_perGene(j);
        end
        if additionalRxns_perGene(j) > 0
            additionalRxns = setdiff(relatedRxns_i, relatedRxns_manual);
            eqs_additional = getRxn_cobraFormat(model_i, additionalRxns);
            st_1 = strcat(additionalRxns, ':' ,eqs_additional);
            st2 = strjoin (st_1,', ');
            st2 = strcat('. Additional: ', st2);
            percentagesAdditional{j,i} = additionalRxns_perGene(j);
        end
        gene_rxns_Matrix{j,i} = strcat(st1, st2);
    end
end

save(['gene_rxns_Matrix_' species], 'gene_rxns_Matrix');
save(['percentagesNotRecovered_' species], 'percentagesNotRecovered');
save(['percentagesAdditional_' species], 'percentagesAdditional');

xlswrite(['gene_rxns_Matrix_' species '.xlsx'],gene_rxns_Matrix)
end