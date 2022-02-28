function createInputsForcEMAF(model, biomassRxnID, baseDir, modelName, constraints_ids, constraints_lb, constraints_ub, option)

changeCobraSolver('gurobi')
cd(baseDir);
fidg = fopen(['EMAF_' modelName '.bat'],'w+');

%it is assumed that reactions in the backward sense will be in the form '
%-> 1 fru[e]'. For that all the exchange reactions in the original model
%should be in the form "1 fru[e] -> "
posSharedExchangeUptake = findExcRxnsWithIDs(model, option);
EXiRxns = model.rxns(posSharedExchangeUptake);
if option==4
    EXiRxns = union(EXiRxns, model.rxns(find(~cellfun(@isempty, regexp(model.rxns,'^spe._')))));
end

if exist([baseDir filesep 'emaf'], 'dir')~=7; mkdir([baseDir filesep 'emaf']); end
if exist([baseDir filesep 'emaf' filesep 'lpx'], 'dir')~=7; mkdir([baseDir filesep 'emaf' filesep 'lpx']); end
if exist([baseDir filesep 'emaf' filesep 'models'], 'dir')~=7; mkdir([baseDir filesep 'emaf' filesep 'models']); end
if exist([baseDir filesep 'emaf' filesep 'minimizationsets'], 'dir')~=7; mkdir([baseDir filesep 'emaf' filesep 'minimizationsets']); end
if exist([baseDir filesep 'emaf' filesep 'constraints'], 'dir')~=7; mkdir([baseDir filesep 'emaf' filesep 'constraints']); end

cd([baseDir filesep 'emaf' filesep 'models'])
newName = [modelName '_irrev.xml'];
if exist(newName,'file')~=2
    writeCbModel(model, 'fileName', 'model', 'format', 'sbml');
    removeSquareBracketsFromBIGGModel('model.xml', newName);
    delete('model.xml');
end

cd([baseDir filesep 'emaf' filesep 'minimizationsets'])
if isempty(find(~cellfun(@isempty, regexp(EXiRxns, '^R_'))))
    EXiRxns = strcat('R_', EXiRxns);
end
data = ['Reaction ID'; EXiRxns];
writeCellArrayInCSVFormat(data, 'all.csv');

cd([baseDir filesep 'emaf' filesep 'constraints'])
if isempty(find(~cellfun(@isempty, regexp(constraints_ids, '^R_'))))
    constraints_ids = strcat('R_', constraints_ids);
end
info = [constraints_ids, num2cell(constraints_lb), num2cell(constraints_ub)];
for i = 1:size(info,1);  for j = 2:size(info,2); info{i,j} = num2str(info{i,j}); end; end;
data = [[{'Reaction ID'},{'LB'},{'UB'}]; info];
writeCellArrayInCSVFormat(data, 'basic.csv');

setModelNameForEMAF(newName, [baseDir filesep 'emaf']);

filePath = regexprep([baseDir filesep 'emaf'],'\\','\\\\');
fprintf(fidg, ['cd /d "' filePath '"\n']);
filePath = regexprep([baseDir filesep 'emaf' filesep 'runMedia3.py'],'\\','\\\\');
fprintf(fidg, ['python "' filePath '"\n']);
filePath = regexprep([baseDir filesep 'emaf' filesep 'pushRunMedia3.py'],'\\','\\\\');
fprintf(fidg, ['python "' filePath '"\n']);

fclose(fidg);

end
