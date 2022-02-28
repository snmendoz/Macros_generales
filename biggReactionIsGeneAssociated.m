function is = biggReactionIsGeneAssociated(rule)

rules_in_model = split(rule,'/');
hasRuleInModel = zeros(size(rules_in_model));
for i = 1:length(rules_in_model)
    pos_colon = strfind(rules_in_model{i},':');
    pos_quotation_marks = strfind(rules_in_model{i},'"');
    pos_quotation_marks = pos_quotation_marks(find(pos_quotation_marks>pos_colon));
    rule = rules_in_model{i}(pos_quotation_marks(1)+1:pos_quotation_marks(2)-1);
    if ~isempty(rule)
        hasRuleInModel(i)= 1;
    end
end

is = any(hasRuleInModel);

end