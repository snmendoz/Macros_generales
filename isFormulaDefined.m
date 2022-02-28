function bool = isFormulaDefined(formula)

bool = 1;
if ~isempty(strfind(formula,'R')) || ~isempty(strfind(formula,'X')) 
   bool = 0; 
end

end