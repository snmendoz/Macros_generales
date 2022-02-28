function printRxnEquationsFromMet(model, mets,names)

if nargin < 3
   names = 0; 
end
 
if ischar(mets)
    pos = find(strcmp(model.mets, mets));
elseif iscell(mets)
    pos = cell2mat(arrayfun(@(x)find(strcmp(x,model.mets)),mets,'UniformOutput',false))';
else
    pos = mets;
end

n_mets = length(pos);
listRxns = [];

if iscell(mets)
    equationsByMet = cell(size(mets));
    for i = 1:n_mets
        posRxns = find(model.S(pos(i), :));
        equationsByMet{i} = getRxn_cobraFormat(model, posRxns,names);
        listRxns = union(listRxns, posRxns);
    end
else
    posRxns = find(model.S(pos, :));
    allEquations = getRxn_cobraFormat(model, posRxns,names);
end

allEquations = strcat(model.rxns(posRxns),' : ',allEquations);
for i = 1:size(allEquations,1)
    fprintf([allEquations{i} '\n']);
end
end