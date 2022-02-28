function exportNotUniformArrayToExcel(info,model,excelFileName)

table = cell(5000,4);
posStart = 1;
for i = 1:length(info)
    if posStart > 1
        posEmpty = find(cellfun(@isempty,table(:,3)));
        posStart = posEmpty(1);
    end
    rxns = info{i,2};
    pos_rxns = getPosOfElementsInArray(rxns, model.rxns);
    eqs = getRxn_cobraFormat(model, pos_rxns);
    rules = model.grRules(pos_rxns);
    n_rxns = length(rxns);
    posEnd = posStart + n_rxns-1;
    table(posStart,1) = info(i,1);    
    table(posStart:posEnd,2) = rxns;
    table(posStart:posEnd,3) = eqs;
    table(posStart:posEnd,4) = rules;
    posStart = posEnd+1;
end
table = table(1:posEnd,:);
xlswrite( excelFileName,table)
end