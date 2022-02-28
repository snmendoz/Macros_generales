function satisfyRule = verifyRule(rule, genes)

genes_Involucrados_por_Rxn=regexp(regexprep(rule,'or|and|\(|\)',' '),'\ ','split');
genes_Involucrados_por_Rxn = setdiff(genes_Involucrados_por_Rxn, '');

boleanRule = rule;
for i = 1:length(genes) 
    boleanRule = regexprep(boleanRule,['(' genes{i} ' '],'1');
    boleanRule = regexprep(boleanRule,[' ' genes{i} ' '],'1');
    boleanRule = regexprep(boleanRule,[' ' genes{i} ')'],'1');
    if isempty(strfind(boleanRule,' ')) && isempty(strfind(boleanRule,'(')) && isempty(strfind(boleanRule,')'))
        boleanRule = regexprep(boleanRule,genes{i},'1');
    end
end

diff = setdiff(genes_Involucrados_por_Rxn, genes);
for i = 1:length(diff) 
    boleanRule = regexprep(boleanRule,['(' diff{i} ' '],'0');
    boleanRule = regexprep(boleanRule,[' ' diff{i} ' '],'0');
    boleanRule = regexprep(boleanRule,[' ' diff{i} ')'],'0');
    if isempty(strfind(boleanRule,' ')) && isempty(strfind(boleanRule,'(')) && isempty(strfind(boleanRule,')'))
        boleanRule = regexprep(boleanRule,diff{i},'0');
    end
end

boleanRule = regexprep(boleanRule, 'and','&');
boleanRule = regexprep(boleanRule, 'or','|');
satisfyRule = 0;
if eval(boleanRule)
    satisfyRule = 1;
end

end