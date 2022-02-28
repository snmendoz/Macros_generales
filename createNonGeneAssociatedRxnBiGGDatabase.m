function createNonGeneAssociatedRxnBiGGDatabase

load('D:\Dropbox\Databases\BIGG\models.mat')
load('D:\Dropbox\Databases\BIGG\bigg_85.mat')
load('D:\Dropbox\Databases\BIGG\Presence_rxns_bigg_85.mat')

% for i = 1:length(models)
%     if ~isfield(models{i}, 'rules') && isfield(models{i}, 'grRules')
%        models{i} = createRulesFromgrRules(models{i}); 
%     end
% end
% 
% withoutGeneAssociations = find(~cellfun(@(x) isfield(x, 'rules'), models));

%let create a matrix of n(rxns) x m(models). 1 if reaction i is has any
%gene associated in model j
hasGeneAssociation = zeros(size(Presence_rxns));
for i = 1:length(bigg.rxns)
    fprintf(['i: ' num2str(i) '\n']);
    inModels = find(Presence_rxns(i,:));
    for j = 1:length(inModels)
        rule_i_j = models{inModels(j)}.grRules{strcmp(models{inModels(j)}.rxns, bigg.rxns{i})==1};
        if ~isempty(rule_i_j)
            hasGeneAssociation(i,inModels(j)) = 1;
        end    
    end
end

rxnsWithoutGeneAssociationsInAllModels = find(sum(hasGeneAssociation,2)==0);

end