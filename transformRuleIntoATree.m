function tree = transformRuleIntoATree(rule)

rule = regexprep(rule, ' ', '');
rule = removeUnnecessaryParenthesisInRule(rule);
posOpen = strfind(rule,'(');
posClose = strfind(rule,')');
level_o = zeros(length(posOpen),1);
level_c = zeros(length(posOpen),1);
for i = 1:length(posOpen)
    if i ==1; level_o(i)=1;
    else
        level_o(i) = length(find(posOpen<posOpen(i))) - length(find(posClose<posOpen(i))) + 1;
    end
end

for i = length(posClose):-1:1
    if i == length(posClose); level_c(i) = 1;
    else
        level_c(i) = length(find(posOpen<posClose(i))) - length(find(posClose<=posClose(i))) + 1;
    end
end
pairs = zeros(length(posOpen),2);
for i = 1:length(posOpen) 
    pairs(i,1) = posOpen(i);
    pos = intersect(find(level_c == level_o(i)),find(posClose>posOpen(i)));
    pos = pos(1);
    pairs(i,2) = posClose(pos);
end

pos_or = strfind(rule, 'or');
pos_and = strfind(rule, 'and');
level_or = zeros(size(pos_or));
level_and = zeros(size(pos_and));
for i = 1:length(pos_or)
    level_or(i) = length(find(posOpen<pos_or(i))) - length(find(posClose<pos_or(i))) + 1;
end
for i = 1:length(pos_and)
    level_and(i) = length(find(posOpen<pos_and(i))) - length(find(posClose<pos_and(i))) + 1;
end

for i = 1:length(pos_or)
    rule = [rule(1:pos_or(i)-1) '|' num2str(level_or(i)) rule(pos_or(i)+2:end)];
end

for i = 1:length(pos_and)
    rule = [rule(1:pos_and(i)-1) '&&' num2str(level_and(i)) rule(pos_and(i)+3:end)];
end
rule = regexprep(rule,'&&','&');

tree = getSetFromRule(rule);

end