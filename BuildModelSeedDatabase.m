function BuildModelSeedDatabase

[n,s] = xlsread('D:\Dropbox\Databases\ModelSEEDDatabase-master\Biochemistry\reactions.xlsx');
[n2,s2] = xlsread('D:\Dropbox\Databases\ModelSEEDDatabase-master\Biochemistry\compounds.xlsx');

load('D:\Projects\BIGG\e_coli_core.mat')
model = e_coli_core;
model.rxns = strcat('OWN_', model.rxns);
rxns = model.rxns;

metaboliteIDs = s2(2:end,1);
metaboliteNames = s2(2:end,3);
metaboliteFormula = s2(2:end,4);
metaboliteCharge = n2(1:end,4);

newIDs = s(2:end,1);
rxnNames = s(2:end,3);
formulas = s(2:end,7); 
formulas = regexprep(formulas, {'(',')','<=','=>','<->','\[0\]','\[1\]'},{'','','<-','->','<=>','[c]','[e]'});
skippedRxns = find(cellfun(@isempty,strfind(formulas,'[2]'))==0);

for i = 1:length(s)-1
    disp(i)
    if ismember(i,skippedRxns); continue; end;
    
    if strfind(formulas{i}, '<-'); 
        pos = strfind(formulas{i}, '<-');
        newFormula = [formulas{i}(pos+3:end) ' -> ' formulas{i}(1:pos-2)];
        formulas{i} = newFormula;
    end
    model = addReaction(model,newIDs{i},'reactionFormula',formulas{i});
end

posa = cell2mat(arrayfun(@(x)find(strcmp(x,newIDs)),model.rxns,'UniformOutput',false))';
model.rxnNames = rxnNames(posa);

model = removeRxns(model, rxns);
mets = model.mets;
mets = regexprep(mets,{'\[c\]','\[e\]'},{'',''});
pos = cell2mat(arrayfun(@(x)find(strcmp(x,metaboliteIDs)),mets,'UniformOutput',false))';

model.metNames = metaboliteNames(pos);
model.metFormulas = metaboliteFormula(pos);
model.metCharge = metaboliteCharge(pos);
model = rmfield(model, 'genes');
model = rmfield(model, 'rules');
model = rmfield(model, 'rev');
model = rmfield(model, 'grRules');
model = rmfield(model, 'rxnGeneMat');
model.csense = repmat('E',length(model.mets),1);
model.description = 'ModelSEED';

save('ModelSEED','model');
writeCbModel(model, 'format', 'sbml', 'fileName', 'ModelSEED');

end