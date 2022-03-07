open function tf = startsWith(str, suffix)

tf = false;
if ~isempty(regexp(str,['^' suffix]))
    tf = true;
end

end