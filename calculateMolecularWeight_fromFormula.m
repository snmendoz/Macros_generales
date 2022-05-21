function mw = calculateMolecularWeight_fromFormula(metFormula)

% List of Elements
elements = {'H','C', 'O', 'P', 'S', 'N', 'Mg', 'Mn', 'Mo', 'Fe','Zn','Co','Ca','I',...
    'Na','Cl','K','Ag','As','Cd','Cu','Hg','Ni','Se','X','R','FULLR'};
molecularWeigths = [1.00784, 12.0107, 15.999, 30.973762, 32.065, 14.0067, ...
    24.305, 54.938044, 95.95, 55.845, 65.38, 58.933195, 40.078,  126.90447, ...
    22.989769, 35.453, 39.0983, 107.8682, 74.9216, 112.411, 63.546,...
    200.59, 58.6934, 78.96 ];
E = zeros(length(elements),1);

% metFormula = model.metFormulas{getPosOfElementsInArray({met}, model.mets)};

for i = 1:length(elements)
    E(i) = numAtomsOfElementInFormula(metFormula, elements{i}, 0);
end
if any(E(25:27)>0)
    mw = 0;
else
    mw = sum(E(1:24).*molecularWeigths');
end

end