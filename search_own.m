function [match ,pos] = search_own(str, strArray)

pos = find(~cellfun(@isempty, strfind(strArray, str)));
match = strArray(pos);

end