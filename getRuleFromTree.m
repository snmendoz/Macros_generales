function rule = getRuleFromTree(tree)

if ischar(tree)
    rule = tree;
else
    if length(tree)==1
        if ischar(tree{1})
            rule = tree{1};
        elseif iscell(tree{1}) && length(tree{1}) == 1
            rule = tree{1};
        elseif iscell(tree{1}) && length(tree{1}) > 1
            rule = strjoin(tree{1}, ' & ');
        end
    else
        strings = cellfun(@getRuleFromTree, tree, 'UniformOutput', 0);
        rule = strjoin(strings, ' | ');
    end
end


end