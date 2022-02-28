function [allMetabolites, metabolites_perEquation] = getMetabolitesNotInModel_fromEquations(model, equations, relaxCompartments)

if nargin<3
    relaxCompartments = 0;
end

allMetabolites = {};
metabolites_perEquation = {};

if ischar(equations)
    metaboliteList = parseRxnFormula(equations);
    if relaxCompartments
        diff = setdiff(removeCompartmentFromMets(metaboliteList),removeCompartmentFromMets(model.mets));
    else
        diff = setdiff(metaboliteList,model.mets);
    end
    if ~isempty(diff)
        metabolites_perEquation = diff;
        allMetabolites = diff;
    end
    
elseif iscell(equations)
    metabolites_perEquation = cell(size(equations));
    for i =1:length(equations)
        metaboliteList = parseRxnFormula(equations{i});
        if relaxCompartments
            diff = setdiff(removeCompartmentFromMets(metaboliteList),removeCompartmentFromMets(model.mets));
        else
            diff = setdiff(metaboliteList,model.mets);
        end
        if ~isempty(diff)
            metabolites_perEquation{i} = diff;
            allMetabolites = union(allMetabolites,diff);
        end
    end
end

end