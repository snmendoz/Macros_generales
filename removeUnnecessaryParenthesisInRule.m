function [newRule, posOpen, posClose, level_o, level_c, pairs] = removeUnnecessaryParenthesisInRule(rule)

newRule = rule;
%remove unnecessary parenthesis
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

if length(posOpen) == length(posClose)
    if ismember(1,posOpen) && ismember(length(rule),posClose)
        if pairs(1,1)==1 && pairs(1,2)==length(rule)
           newRule = rule(2:end-1);
        end
    end
else
    error('wrong parenthesis')
end

end