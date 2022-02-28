function mw = calculateMolecularWeight(model,met)

% List of Elements
elements = {'H','C', 'O', 'P', 'S', 'N', 'Mg','Fe','Zn','Co','Ca','Y','I','Na','Cl','K','X','R','FULLR'};
molecularWeigths = [1.00784, 12.0107, 15.999, 30.973762, 32.065, 14.0067, 24.305, 55.845, 65.38, 58.933195, 40.078, 88.90585, 126.90447, 22.989769, 35.453, 39.0983];
E = zeros(length(elements),1);

metFormula = model.metFormulas{getPosOfElementsInArray({met}, model.mets)};

for i = 1:length(elements)
    E(i) = numAtomsOfElementInFormula(metFormula, elements{i}, 0);
end
if any(E(17:19)>0)
    mw = 0;
else
    mw = sum(E(1:16).*molecularWeigths');
end

end