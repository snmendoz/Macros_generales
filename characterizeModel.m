function characterizeModel(model)
exchange = model.rxns(findExcRxns(model));
transport = findTransRxns(transformModelToCOBRAFormat( model));
cytosol = setdiff(model.rxns,union(transport,exchange));

gene_associated  = model.rxns(find(~cellfun(@isempty, model.grRules)));
gene_associated_cytosol = intersect(gene_associated,cytosol);
gene_associated_transport = intersect(gene_associated,transport);

not_gene_associated = setdiff(model.rxns,gene_associated);
not_gene_associated_cytosol = intersect(not_gene_associated,cytosol);
not_gene_associated_transport = intersect(not_gene_associated,transport);
not_gene_associated_exchange = intersect(not_gene_associated,exchange);

labels = [[{'gene_associated'};{''};{''};{'not_gene_associated'};{''};{''};{''}],...
    [{''};{'gene_associated_cytosol'};{'gene_associated_transport'};{''};{'not_gene_associated_cytosol'};{'not_gene_associated_transport'};{'not_gene_associated_exchange'}]];
numbers = [length(gene_associated);length(gene_associated_cytosol);length(gene_associated_transport);length(not_gene_associated);...
    length(not_gene_associated_cytosol);length(not_gene_associated_transport);length(not_gene_associated_exchange)];
info = [labels,num2cell(numbers)];
xlswrite('basic_info',info,'basic_info');


percentage_gene_associated = 100*(length(gene_associated))/(length(gene_associated)+length(not_gene_associated_cytosol)+length(not_gene_associated_transport));
xlswrite('basic_info',percentage_gene_associated,'percentage_ga');

info = [{'n_rxns'};{'n_mets'};{'n_genes'}];
info = [info, num2cell([length(model.rxns);length(model.mets);length(model.genes)])];
xlswrite('basic_info',info,'general');
end