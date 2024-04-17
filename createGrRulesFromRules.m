function model = createGrRulesFromRules(model)

rules = model.rules;
grRules = rules;

involvedGenesPerReaction = cellfun(@unique,cellfun(@splitString,regexprep(grRules,{'\(|\)','\ or\ |\ and\ |\ &\ |\ \|\ '},{'',' '}),'UniformOutput',0),'UniformOutput',0);
% involvedGenesPerReaction=regexp(regexprep(model.rules,'\||&|\(|\)',''),'\  ','split');
% for i = 1:length(involvedGenesPerReaction); 
%     for j = 1:length(involvedGenesPerReaction{i})
%         involvedGenesPerReaction{i}{j} = regexprep(involvedGenesPerReaction{i}{j},'\ ','');
%     end
% end

for i =1:length(model.rxns)
    for j = 1:length(involvedGenesPerReaction{i})
%         fprintf(['i:' num2str(i) ' j:' num2str(j) '\n']);
        if ~isempty(involvedGenesPerReaction{i}{j})
            grRules{i} = regexprep(grRules{i}, [involvedGenesPerReaction{i}{j}(1) '\(' involvedGenesPerReaction{i}{j}(2:end) '\)'], model.genes(str2num(involvedGenesPerReaction{i}{j}(2:end))));
        end
    end
end

grRules = regexprep(grRules,  ' \| ',' or ');
grRules = regexprep(grRules,  ' \& ',' and ');
model.grRules = grRules;
end