function model = setMediaFromRxns(model, mediaExchangeRxns, valueMediaExchangeRxns,allNutrientExchangeRxns)

if all(arrayfun(@(x) full(model.S(find(model.S(:,x)), x))<0, getPosOfElementsInArray(allNutrientExchangeRxns,model.rxns)))
    model = changeRxnBounds(model,allNutrientExchangeRxns,0,'l');
    model = changeRxnBounds(model,mediaExchangeRxns,valueMediaExchangeRxns,'l');
elseif all(arrayfun(@(x) full(model.S(find(model.S(:,x)), x))<0, getPosOfElementsInArray(allNutrientExchangeRxns,model.rxns))==0)
    model = changeRxnBounds(model,allNutrientExchangeRxns,0,'u');
    model = changeRxnBounds(model, mediaExchangeRxns(ismember(mediaExchangeRxns,model.rxns)),...
        valueMediaExchangeRxns(ismember(mediaExchangeRxns,model.rxns)),'u');

end


end