function tf = contains(str, suffix)
tf = false;
if ~isempty(strfind(str,suffix))
    tf = true;
end

end