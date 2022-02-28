function set = getSetFromRule(rule)
rule = removeUnnecessaryParenthesisInRule(rule);
pos_or = strfind(rule,'|');
pos_and = strfind(rule,'&');

if isempty(pos_or) && isempty(pos_and)
    set = {rule};
elseif ~isempty(pos_or) && isempty(pos_and)
    levels = str2num(rule(pos_or+1)');
    if length(levels) == 1
        set = regexp(rule,['\|' num2str(levels)],'split');
    else
        minLevel = min(levels);
        groups = regexp(rule,['\|' num2str(minLevel)],'split');
        set = cell(length(groups),1);
        for i = 1:length(groups)
            group_i = getSetFromRule(groups{i});
            set{i} = group_i{1};
        end
    end
elseif ~isempty(pos_and) && isempty(pos_or)
    levels = unique(str2num(rule(pos_and+1)'));
    if length(levels) ==1
    set = cell(1,1);
    genes = regexp(rule,['\&' num2str(levels)],'split');
    set{1} = genes;
    else
        error('wrong')
    end
elseif ~isempty(pos_and) && ~isempty(pos_or)
    minLevel_or = min(str2num(rule(pos_or+1)'));
    minLevel_and = min(str2num(rule(pos_and+1)'));
    
    if minLevel_or < minLevel_and
        groups = regexp(rule,['\|' num2str(minLevel_or)],'split');
        set = cell(length(groups),1);
        for i = 1:length(groups)
            if ~isempty(strfind(groups{i},'|')) || ~isempty(strfind(groups{i},'&'))
                group_i = getSetFromRule(groups{i});
            else
                group_i = groups(i);
            end
            set(i) = group_i;
        end
    else
        set = cell(1,1);
        genes = regexp(rule,['\&' num2str(minLevel_and)],'split');
        set{1} = genes;
        for j = 1:length(genes)
            if ~isempty(strfind(genes{j},'|')) || ~isempty(strfind(genes{j},'&'))
                set{1}{j} = getSetFromRule(genes{j});
            end
        end
    end
    
end

end