function exportMassBalanceToCSV(model, fileName, includeBoundaryMetabolites)


info = cell(length(model.rxns),2);
info(:,1) = model.rxns;
excRxns = find(findExcRxns(model));
for i = 1:length(model.rxns)
    balance = '';
    stoi = model.S(find(model.S(:,i)),i);
    mets = model.mets(find(model.S(:,i)));
    
    reactants_ids = mets(find(stoi<0)); 
    reactants_stoi = num2cell(abs(stoi(find(stoi<0))));
    for j = 1:length(reactants_stoi)
        reactants_stoi{j} = num2str(reactants_stoi{j});
    end
    
    if includeBoundaryMetabolites && ismember(i, excRxns)
        products_ids = transformMetsToCBMPYFormat(regexprep(transformMetsToCOBRAFormat(reactants_ids),'\[.*\]','[boundary]'));
        products_stoi = {'1'};
    else
        products_ids = mets(find(stoi>0));
        products_stoi = num2cell(stoi(find(stoi>0)));
        for j = 1:length(products_stoi)
            products_stoi{j} = num2str(products_stoi{j});
        end
        
    end
    for j = 1:length(reactants_stoi)
       balance = [ balance '- ' reactants_stoi{j} ' ' reactants_ids{j} ' ']; 
    end
    for j = 1:length(products_stoi)
       balance = [ balance '+ ' products_stoi{j} ' ' products_ids{j} ' ']; 
    end
    balance = balance(1:end-1);
    
    info{i,2} = balance;
end
exportToCSV(fileName, info)

end