function [rxnFormulaWasFound, posRxn, rxnOppositeDirection, matchWithoutRemovingProtons, additionalMatchsIfRemovingProtons, posRxnsFullMatch] = reactionFormulaInModel2(model, equation, allowProductsByReactants, allowProtonMismatch)

if nargin < 4
    allowProtonMismatch = 0;
end

rxnFormulaWasFound = 0;
posRxn = [];
posRxnsFullMatch = [];
rxnOppositeDirection = 0;
matchWithoutRemovingProtons = 0;
additionalMatchsIfRemovingProtons = 0;

if ischar(equation)
    [metaboliteList, stoichCoeffList, ~] = parseRxnFormula(equation);
elseif iscell(equation)
    [metaboliteList, stoichCoeffList, ~] = parseRxnFormula(equation{1});
elseif isnumeric(equation)
    posMets = find(model.S(:,equation));
    metaboliteList = model.mets(posMets);
    stoichCoeffList = full(model.S(posMets, equation));
elseif isstruct(equation)
    modelRef = equation.model;
    posRxnRef = find(strcmp(modelRef.rxns,equation.rxn));
    posMets = find(modelRef.S(:,posRxnRef));
    metaboliteList = modelRef.mets(posMets);
    stoichCoeffList = full(modelRef.S(posMets, posRxnRef));
end

[isInModel, posMets] = ismember(metaboliteList, model.mets);

if any(~isInModel)
    return
else
    
    int = find(model.S(posMets(1),:));
    intersectedMets = 1;
    n_restantes = length(posMets) - intersectedMets;
    
    while n_restantes > 0 && ~isempty(int)
        intersectedMets = intersectedMets + 1;
        n_restantes = length(posMets) - intersectedMets;
        int = intersect(int, find(model.S(posMets(intersectedMets),:)));
    end
    if ~isempty(int)
        int = int(arrayfun(@(x) length(find(model.S(:,x))), int) == length(posMets));
        posRxn = zeros(length(int),1);
        rxnOppositeDirection = zeros(length(int),1);
        cont = 0;
        
        for i = 1:length(int)
            stoi = full(model.S(posMets,int(i)));
            if size(stoichCoeffList,1) ~= size(stoi,1); stoichCoeffList = stoichCoeffList'; end;
            if all(stoi == stoichCoeffList)
                cont = cont+1;
                posRxn(cont) = int(i);
            elseif (all(stoi == -stoichCoeffList) && allowProductsByReactants)
                cont = cont+1;
                rxnOppositeDirection(cont) = 1;
                posRxn(cont) = int(i);
            end
        end
        
        if cont > 0
            rxnFormulaWasFound = 1;
            posRxn = posRxn(1:cont);
            rxnOppositeDirection = rxnOppositeDirection(1:cont);
        end
    end
    
    if alloallowProtonMismatch
        isProtonInFormula = ~isempty(find(~cellfun(@isempty, regexp(metaboliteList, '^h_.*$'))));
        if isProtonInFormula
            posProtons = find(~cellfun(@isempty, regexp(metaboliteList, '^h_.*$')));
            
            %first for reactions with different proton stoichiometry 
            int = find(model.S(posMets(1),:));
            intersectedMets = 1;
            n_restantes = length(posMets) - intersectedMets;
            
            while n_restantes > 0 && ~isempty(int)
                intersectedMets = intersectedMets + 1;
                n_restantes = length(posMets) - intersectedMets;
                int = intersect(int, find(model.S(posMets(intersectedMets),:)));
            end
            if ~isempty(int)
                int = int(arrayfun(@(x) length(find(model.S(:,x))), int) == length(posMets));
                posRxn = zeros(length(int),1);
                rxnOppositeDirection = zeros(length(int),1);
                cont = 0;
                
                for i = 1:length(int)
                    stoi = full(model.S(posMets,int(i)));
                    if size(stoichCoeffList,1) ~= size(stoi,1); stoichCoeffList = stoichCoeffList'; end;
                    
                    allEqualWithAnalazingProtons = 1;
                    for j = 1:length(stoi)
                    end
%                     if all(stoi == stoichCoeffList)
%                         cont = cont+1;
%                         posRxn(cont) = int(i);
%                     elseif (all(stoi == -stoichCoeffList) && allowProductsByReactants)
%                         cont = cont+1;
%                         rxnOppositeDirection(cont) = 1;
%                         posRxn(cont) = int(i);
%                     end
                end
                
                if cont > 0
                    rxnFormulaWasFound = 1;
                    posRxn = posRxn(1:cont);
                    rxnOppositeDirection = rxnOppositeDirection(1:cont);
                end
            end
            
            %second, for reactions with no protons
            
            kept = true(size(metaboliteList)); kept(posProtons) = false;
            posMets = posMets(kept);
            stoichCoeffList = stoichCoeffList(kept);
            metaboliteList = model.mets(posMets);
            posProtonsInModel = find(~cellfun(@isempty, regexp(model.mets, '^h_.*$')));
            model = removeMetabolites(model, model.mets(posProtonsInModel));
            posMets = cell2mat(arrayfun(@(x)find(strcmp(x,model.mets)),metaboliteList,'UniformOutput',false))';
        else
        end
    end
    
end


end
