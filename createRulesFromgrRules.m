function model = createRulesFromgrRules(model, createGenesIfNeeded)

if nargin<2
    createGenesIfNeeded = 0;
end

grRules = model.grRules;

for i = 1:length(grRules)
    if isempty(grRules{i})
        grRules{i} = '';
    end
end
rules = grRules;
involvedGenes = cellfun(@unique,cellfun(@splitString,regexprep(grRules,{'\(|\)','\ or\ |\ and\ |\ &\ |\ \|\ '},{'',' '}),'UniformOutput',0),'UniformOutput',0);

if createGenesIfNeeded
    allGenes = {};
    for i = 1:length(involvedGenes)
        allGenes = union(allGenes,involvedGenes{i});
    end
end
% involvedGenes = cellfun(@unique,cellfun(@splitString,regexprep(grRules,'\or|and|\(|\)|&|\|',' '),'UniformOutput',0),'UniformOutput',0);
pos = find(~cellfun(@isempty, involvedGenes));
for i = 1:length(pos)
%     if pos(i) ==410
%        disp('') 
%     end
    genes_i = unique(involvedGenes{pos(i)});
%     if pos(i)==getPosOfElementsInArray({'TRE6PP'},model.rxns)
%         disp('')
%     end
    for j = 1:length(genes_i)
        condition1 = ['^' genes_i{j} '$' ... %YALI0F04$
        '|^(' genes_i{j} ')$'...%^(YALI0F04)$    
        '|^' genes_i{j} '\ '... %^YALI0F04\E 
        '|^(' genes_i{j} '\ '....%^(YALI0F04\E 
        '|^(' genes_i{j} ')\ ' ... %^(YALI0F04)E
        '|\ ' genes_i{j} '$'... %\ YALI0F04$
        '|\ (' genes_i{j} ')$'...  %(YALI0F04)$
        '|\ ' genes_i{j} ')$'...  %YALI0F04)$
        '|\ ' genes_i{j} '\ '... %\ YALI0F04\E 
        '|^\ (' genes_i{j} '\ '....%\E(YALI0F04\E
        '|\ ' genes_i{j} ')\ '... %\EYALI0F04)\E 
        '|\ (' genes_i{j} ')\ ']; %\E(YALI0F04)\E 
    
        conditions = {['^' genes_i{j} '$'] ... %YALI0F04$
        ['^\(' genes_i{j} '\)$']...%^(YALI0F04)$    
        ['^' genes_i{j} '(?=\ )']... %^YALI0F04\E 
        ['^\(' genes_i{j} '(?=\ )']....%^(YALI0F04\E 
        ['^\(' genes_i{j} '\)(?=\ )'] ... %^(YALI0F04)E
        ['(?<=\ )' genes_i{j} '$']... %\ YALI0F04$
        ['(?<=\ )\(' genes_i{j} '\)$']...  %(YALI0F04)$
        ['(?<=\ )' genes_i{j} '\)$']...  %YALI0F04)$
        ['(?<=\ )' genes_i{j} '(?=\ )']... %\ YALI0F04\E 
        ['(?<=\ )\(' genes_i{j} '(?=\ )']....%\E(YALI0F04\E
        ['(?<=\ )' genes_i{j} '\)(?=\ )']... %\EYALI0F04)\E 
        ['(?<=\ )\(' genes_i{j} '\)(?=\ )']...
        ['(?<=\()\(' genes_i{j} '\)(?=\ )']...
        ['(?<=\ )\(' genes_i{j} '\)(?=\))']...
        ['(?<=\()\(' genes_i{j} '\)(?=\))']...
        ['^\(\(' genes_i{j} '(?=\ )'],...
        ['^\(\(\(' genes_i{j} '(?=\ )'],...
        ['(?<=\ )' genes_i{j} '\)\)$'],...
        ['(?<=\ )' genes_i{j} '\)\)\)$'],...
        ['\(\(' genes_i{j} '(?=\ )'],...
        ['(?<=\ )' genes_i{j} '\)\)']}; 
%       
        for k =1:length(conditions)
%             if pos(i) ==410
%                 regexprep(rules{pos(i)}, conditions{k} , ['x(' num2str(find(strcmp(genes_i{j}, model.genes))) ')'])
%             end
            rules{pos(i)} = regexprep(rules{pos(i)}, conditions{k} , ['x(' num2str(find(strcmp(genes_i{j}, model.genes))) ')']);
        end
%         condition = ['^' genes_i{j} '$|^' genes_i{j} '(?=\ )|(?<=\ )' genes_i{j} '$|(?<=\ )'  genes_i{j} '(?=\ )|(' genes_i{j} '(?=\ )|(?<=\ )' genes_i{j} ')|(' genes_i{j} ')'];
%         rules{pos(i)} = regexprep(rules{pos(i)}, condition , ['x(' num2str(find(strcmp(genes_i{j}, model.genes))) ')']);
    end
end

rules = regexprep(rules, ' or ', ' | ');
rules = regexprep(rules, ' and ', ' & ');
model.rules = rules;

end