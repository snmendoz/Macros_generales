function is = isMetaboliteInRxnEquation(met,eq)

if iscell(eq); eq = eq{1}; end
is = ismember(met, parseRxnFormula(eq));

end