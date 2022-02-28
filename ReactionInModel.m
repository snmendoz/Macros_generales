function [rxnInModel, rxnIDexists, rxnOppositeDirection] = reactionFormulaInModel(model, equation, allowProductsByReactants)
rxnIDexists = -1;
rxnInModel = false;
rxnOppositeDirection = 0;

if nargin <3; allowProductsByReactants = 0; end; 

nMets = length(model.mets);
Scolumn = sparse(nMets,1);
[metaboliteList, stoichCoeffList, ~] = parseRxnFormula(equation);
[isInModel, metID] = ismember(metaboliteList, model.mets);
nNewMets = sum(~isInModel);

if nNewMets>0; return; end;

for i = 1:length(metaboliteList)
    Scolumn(metID(i),1) = stoichCoeffList(i);
end

rxnInModel=false;

Stmp = model.S;
if size(Stmp,2) < 6000
    tmpSel = all(repmat((Scolumn),1,size(Stmp,2)) == (Stmp));
    rxnIDexists = full(find(tmpSel));
    if ~isempty(rxnIDexists)
        rxnIDexists = rxnIDexists(1);
        rxnInModel = true;
    end
else
    for i=1:size(Stmp,2)
        if Scolumn == Stmp(:,i)
            rxnInModel = true;
            rxnIDexists = i;
            break
        elseif (-Scolumn == Stmp(:,i))
            if allowProductsByReactants
                rxnOppositeDirection = 1;
                rxnInModel = true;
                rxnIDexists = i;
                break
            end
        end
    end
end

end