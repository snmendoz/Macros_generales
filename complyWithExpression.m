function comply = complyWithExpression(string, expression)

if isempty(string)
    comply = -0;
elseif isempty(cell2mat(regexp(strsplit(string,';'), expression)))
   comply = -1;
elseif all(cell2mat(regexp(strsplit(string,';'), expression)))
    comply = 1;
end

end