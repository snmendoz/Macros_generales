function equal = sameFormula(formula1, formula2)

equal = 0;

elements = {'H','C','N','O','P','S','Na','Mg','Cl','K','Ca','Mn','Fe','Ni','Co','Cu','Zn','As','Se','Ag','Cd','W','Hg','Mo','I','R','X'};
n1 = zeros(length(elements));
n2 = zeros(length(elements));
for i = 1:length(elements)
    n1(i) = numAtomsOfElementInFormula(formula1, elements{i});
    n2(i) = numAtomsOfElementInFormula(formula2, elements{i});
end
if n1==n2;
    equal = 1;
end

end