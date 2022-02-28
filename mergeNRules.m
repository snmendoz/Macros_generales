function rule = mergeNRules(model, pos)

rules = model.grRules(pos);

nRules = length(rules);
if nRules ==2
    rule = mergeRules(rules{1}, rules{2});
else
    rule = mergeRules(rules{1}, rules{2});
    mergedRules = 2;

    rest = nRules-mergedRules;
    while rest > 0
        nextRule = merged + 1;
        rule = mergeRules(rule, rules{nextRule});
        mergedRules = mergedRules + 1;
        rest = nRules - mergedRules;
    end
end

end