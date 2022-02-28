function newRule = removeUnnecessarySpacesInRule(rule)

newRule = rule;
while strcmp(rule(1),' ')
    rule = rule(2:end);    
end

while strcmp(rule(end), ' ')
    rule = rule(1:end-1);    
end

posOpen = strfind(rule,'(');
posClose = strfind(rule,')');
posSpaces = strfind(rule,' ');

while any(ismember(posSpaces,posOpen+1))
    pos = find(ismember(posSpaces,posOpen+1));
    for i = 1:length(pos)
        rule = [rule(1:posSpaces(pos(i))-1) rule(posSpaces(pos(i))+1:end)];
    end
    posOpen = strfind(rule,'(');
    posClose = strfind(rule,')');
    posSpaces = strfind(rule,' ');
end

while ismember(posSpaces,posClose-11)
    posOpen = strfind(rule,'(');
    posClose = strfind(rule,')');
    posSpaces = strfind(rule,' ');
end

end