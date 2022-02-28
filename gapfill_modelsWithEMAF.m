function gapfill_modelsWithEMAF(models, idsType, names, abbr, biomass, metOtherIDs, metMNXIDs, rxnOtherIDs, rxnMNXIDs, species)
fprintf('performing gap-filling\n');
modelRef = models{1};
fbaRef = optimizeCbModel(modelRef);

fidg = fopen('EMAF_gapfilling.bat','w+');

threshold = 10^-6;
if exist('minFlux2.mat','file') ==2 && exist('maxFlux2.mat','file') ==2
    load('minFlux2.mat');
    load('maxFlux2.mat');
else
    minFlux2 = zeros(size(modelRef.rxns));
    maxFlux2 = zeros(size(modelRef.rxns));
    for i = 1:length(modelRef.rxns)
        model2 = changeObjective(modelRef, modelRef.rxns(i));
        fmin = optimizeCbModel(model2, 'min');
        fmax = optimizeCbModel(model2, 'max');
        minFlux2(i) = fmin.f;
        maxFlux2(i) = fmax.f;
    end
    save('minFlux2.mat','minFlux2');
    save('maxFlux2.mat','maxFlux2');
end
posActiveRxns = union(find(minFlux2<-threshold), find(maxFlux2>threshold));
% posBiomass = find(strcmp(modelRef.rxns,regexprep(biomass,'R_','')));
% posActiveRxns = setdiff(posActiveRxns, posBiomass);
activeRxns = modelRef.rxns(posActiveRxns);
infeasibleNetworks = cell(length(models)-1,1);
setsToSearch = cell(length(models)-1,1);

if exist('infeasibleNetworks.mat','file')~=2
    for i = 2:length(models)
        fprintf('preparing model %2.0f ...',(i-1));
        if ~strcmp(names{i},'CarveMe') && ~strcmp(names{i},'ModelSEED')
            if strcmp(idsType{i}, 'bigg')
                model_i = models{i};
            else
                model_i = removeDuplicatedMetabolitesFromModel(models{i});
                model_i = removeEmptyRxns(model_i);
                model_i = translateModelToTargetLanguage(model_i, idsType{i}, 'bigg', models{1}, metOtherIDs, metMNXIDs, rxnOtherIDs, rxnMNXIDs, 0);
            end
            
            %check that identifiers that intersect have the same formula
            [~, pos1, pos2] = intersect(model_i.rxns, modelRef.rxns);
            [~, sameRxnFormula, ~, ~] = ...
                arrayfun(@(x,y) compareRxns(model_i, x, modelRef, y), pos1, pos2);
            if any(sameRxnFormula==0) 
                disp('')
            end

            [diff, posDiff] = setdiff(model_i.rxns, modelRef.rxns);
            eq = getRxn_cobraFormat(model_i, diff);
            modelAux = modelRef;
            for j = 1:length(diff)
                modelAux = addReaction(modelAux, diff{j}, 'reactionFormula', eq{j},'printLevel',0);
                lb = model_i.lb(posDiff(j));
                ub = model_i.ub(posDiff(j));
                pos = find(strcmp(modelAux.rxns,diff{j}));
                modelAux = changeRxnBounds(modelAux, modelAux.rxns(pos), 'l', lb);
                modelAux = changeRxnBounds(modelAux, modelAux.rxns(pos), 'u', ub);
            end
            modelAux = changeObjective(modelAux, regexprep(biomass,'R_',''));
            fba_check = optimizeCbModel(modelAux,'max');
            
            %         rxns_i = model_i.rxns;
            %         diff = setdiff(activeRxns, rxns_i);
            %         eq = getRxn_cobraFormat(modelRef, diff);
            %         for j = 1:length(diff)
            %             model_i = addReaction(model_i, diff{j}, 'reactionFormula', eq{j});
            %         end
            %         for j = 1:length(modelRef.rxns)
            %             lb = modelRef.lb(j);
            %             ub = modelRef.ub(j);
            %             pos = find(strcmp(model_i.rxns,modelRef.rxns{j}));
            %             model_i = changeRxnBounds(model_i, model_i.rxns(pos), 'l', lb);
            %             model_i = changeRxnBounds(model_i, model_i.rxns(pos), 'u', ub);
            %         end
            %         model_i = changeObjective(model_i, regexprep(biomass,'R_',''));
            %         fba_check = optimizeCbModel(model_i,'max');
            if fba_check.f>threshold
                setToSearch = setdiff(modelRef.rxns, model_i.rxns);
                
                [modelIrrev, ~, rev2irrev, ~] = convertToIrreversible(modelAux,'sRxns',setToSearch);
                for k = 1:length(rev2irrev); if length(rev2irrev{k})>1; rev2irrev{k} = rev2irrev{k}'; end; end;
                infeasibleNetworks{i-1} = modelIrrev;
                posSetToSearch = cell2mat(arrayfun(@(x)find(strcmp(x,modelAux.rxns)),setToSearch,'UniformOutput',false))';
                setsToSearch{i-1} = strcat('R_', modelIrrev.rxns(cell2mat(rev2irrev(posSetToSearch))));
                
            else
                disp('')
            end
        end
        fprintf('done: progress %2.0f %%\n',100*(i-1)/(length(models)-1));
    end
    save('infeasibleNetworks','infeasibleNetworks')
    save('setsToSearch','setsToSearch')
else
    load('infeasibleNetworks');
    load('setsToSearch');
end

lb = 10^-6;
ub = 1;
baseDir = ['D:\Dropbox\Research_Projects\Review_reconstruction\comparison\' species 'GapFilling_Networks'];

for i = 1:length(infeasibleNetworks)
    fprintf('preparing model %2.0f ...',i);
    if isempty(infeasibleNetworks{i})
        fprintf('done: progress %2.0f %%\n',100*(i)/(length(models)-1));
        continue;
    end
    folder = abbr{i+1};
    
    if exist([baseDir filesep folder filesep 'emaf'], 'dir')~=7; mkdir([baseDir filesep folder filesep 'emaf']); end
    if exist([baseDir filesep folder filesep 'emaf' filesep 'lpx'], 'dir')~=7; mkdir([baseDir filesep folder filesep 'emaf' filesep 'lpx']); end
    if exist([baseDir filesep folder filesep 'emaf' filesep 'models'], 'dir')~=7; mkdir([baseDir filesep folder filesep 'emaf' filesep 'models']); end
    if exist([baseDir filesep folder filesep 'emaf' filesep 'minimizationsets'], 'dir')~=7; mkdir([baseDir filesep folder filesep 'emaf' filesep 'minimizationsets']); end
    if exist([baseDir filesep folder filesep 'emaf' filesep 'constraints'], 'dir')~=7; mkdir([baseDir filesep folder filesep 'emaf' filesep 'constraints']); end
    
    modelIrrev = infeasibleNetworks{i};
    
    cd([baseDir filesep folder filesep 'emaf' filesep 'models'])
    modelFile = [abbr{i+1} '.xml'];
    newName = regexprep(modelFile, '.xml', '_irrev.xml');
    if exist(newName,'file')~=2
        modelIrrev.mets = regexprep(modelIrrev.mets,'_(.(?!_))+$','[$1]');
        exportBIGGModelToSBML(modelIrrev, newName, 1)
%         writeCbModel(modelIrrev, 'fileName', 'model', 'format', 'sbml');
%         removeSquareBracketsFromBIGGModel('model.xml', newName);
%         delete('model.xml');
    end
    
    cd([baseDir filesep folder filesep 'emaf' filesep 'minimizationsets'])
    data = ['Reaction ID'; setsToSearch{i}];
    writeCellArrayInCSVFormat(data, 'all.csv');
    
    cd([baseDir filesep folder filesep 'emaf' filesep 'constraints'])
    data = [[{'Reaction ID'},{'LB'},{'UB'}]; [{biomass}, {num2str(lb,'%1.1f')}, {num2str(ub,'%1.1f')}]];
    writeCellArrayInCSVFormat(data, 'basic.csv');
    
    setModelNameForEMAF(newName, [baseDir filesep folder filesep 'emaf']);
    
    filePath = regexprep([baseDir filesep folder filesep 'emaf'],'\\','\\\\');
    fprintf(fidg, ['cd /d "' filePath '"\n']);
    filePath = regexprep([baseDir filesep folder filesep 'emaf' filesep 'runMedia3.py'],'\\','\\\\');
    fprintf(fidg, ['python "' filePath '"\n']);
    filePath = regexprep([baseDir filesep folder filesep 'emaf' filesep 'pushRunMedia3.py'],'\\','\\\\');
    fprintf(fidg, ['python "' filePath '"\n']);
    
    fprintf('done: progress %2.0f %%\n',100*(i)/(length(models)-1));
end
fclose(fidg);

end