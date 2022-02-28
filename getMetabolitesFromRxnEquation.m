function metaboliteList = getMetabolitesFromRxnEquation(equation)

if ischar(equation)
    metaboliteList = parseRxnFormula(equation);
elseif iscell(equation)
    metaboliteList = parseRxnFormula(equation{1});
end

end