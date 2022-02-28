function [rxnIDs, rxnNames, equations, geneAssociations] = getRxnsFromSimpheny(file)

fi = fopen(file);

tline = fgetl(fi);
tline = fgetl(fi);
rxnIDs = cell(1000,1);
rxnNames = cell(1000,1);
equations = cell(1000,1);
geneAssociations = cell(1000,1);
cont_rxns = 0;

while ischar(tline) 
    fields = regexp(tline, ',','split');
    cont_rxns = cont_rxns + 1;
    rxnIDs{cont_rxns} = regexprep(fields{1},'"','');
    rxnNames{cont_rxns} = regexprep(fields{2},'"','');
    equations{cont_rxns} = regexprep(fields{3},{'"','(',')'},{'','',''});
    if ~isempty(strfind(equations{cont_rxns},':'))
        comp = equations{cont_rxns}(1:3);
        equations{cont_rxns} = equations{cont_rxns}(7:end);
        
        [metaboliteList, stoichCoeffList, revFlag] = parseRxnFormula(equations{cont_rxns});
        metaboliteList = strcat(metaboliteList, comp);
        formula = createRxnFormulaFromProperties(metaboliteList, stoichCoeffList, revFlag);
        equations{cont_rxns} = formula;
    else 
        equations{cont_rxns} = regexprep(equations{cont_rxns},'--','-');
    end
    geneAssociations{cont_rxns} = regexprep(fields{7},'"','');
    tline = fgetl(fi);
end
rxnIDs = rxnIDs(1:cont_rxns);
rxnNames = rxnNames(1:cont_rxns);
equations = equations(1:cont_rxns);
geneAssociations = geneAssociations(1:cont_rxns);

fclose(fi);

end