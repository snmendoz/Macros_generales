function rule = mergeRules(rule1, rule2)
if iscell(rule1) 
    rule1 = rule1{1};
end
if iscell(rule2) 
    rule2 = rule2{1};
end
if isempty(rule1) && isempty(rule2)
    rule = '';
elseif isempty(rule1) && ~isempty(rule2)
    rule = rule2;
elseif isempty(rule2) && ~isempty(rule2)
    rule = rule1;
else
    tree1 = transformRuleIntoATree(rule1);
    tree2 = transformRuleIntoATree(rule2);
    
    tree = [tree1;tree2];
    if all(cellfun(@iscell, tree)==0)
        tree = unique(tree);
    end
    rule = getRuleFromTree(tree);
    
    
    disp('')
end



end