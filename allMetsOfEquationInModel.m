function true = allMetsOfEquationInModel(model, equation, relaxCompartment)

if nargin <3
    relaxCompartment = 0;
end

if ischar(equation)
    [metaboliteList, ~, ~] = parseRxnFormula(equation);
elseif iscell(equation)
    equation = equation{1};
    [metaboliteList, ~, ~] = parseRxnFormula(equation);
elseif isstruct(equation)
    modelRef = equation.model;
    posRxnRef = find(strcmp(modelRef.rxns,equation.rxn));
    posMets = find(modelRef.S(:,posRxnRef));
    metaboliteList = modelRef.mets(posMets);
end

true = 0;
if allMetsInCBMPYFormat(metaboliteList) && allMetsInCOBRAFormat(model.mets)
    metaboliteList = transformMetsToCOBRAFormat(metaboliteList);
elseif allMetsInCBMPYFormat(model.mets) && allMetsInCOBRAFormat(metaboliteList)
    metaboliteList = transformMetsToCBMPYFormat(metaboliteList);
end

if allMetsInModel(model, metaboliteList, relaxCompartment)
   true = 1; 
end
end