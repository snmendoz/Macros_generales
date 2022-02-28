function buildBIGG
% downloadBIGGDatabase
names = createBIGGdatabase;
names_all = names;
for i = 1:length(names)    
    CreatePresenceAndReversibility(names{i});
    addCHEBItoBIGG(names{i})
    addFormulaAndCharge(names{i})
    n = RefineBiggDatabase(names{i});
    names_all = [names, n];
    for j = 1:length(n) 
        m = generateBIGG_cytosol(n{j});
        names_all = [names, m];
    end
end
save('names_all','names_all')
%step 2

load('names_all')
for i = 1:length(names_all)
    exportBIGGDatabaseToSBML(names_all{i})
end

end