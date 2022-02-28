function exportModelPerPathwayToExcelFile(model, excelFile)

posEmpty = find(cellfun(@isempty, model.subSystems));
model.subSystems(posEmpty) = repmat({'Other'},length(posEmpty),1);
uniquePathway = unique(model.subSystems);
uniquePathwayAll = {};
for i = 1:length(uniquePathway)
    pathways = strsplit(uniquePathway{i}, '; ');
    uniquePathwayAll = union(uniquePathwayAll, pathways);
end

posPerPathway = cell(size(uniquePathwayAll));
for i = 1:length(model.rxns)
    if ~isempty(model.subSystems{i})
        pathways = strsplit(model.subSystems{i}, '; ');
        for j = 1:length(pathways)
            posPathway = find(strcmp(uniquePathwayAll,pathways{j} ));
            posPerPathway{posPathway} = union(posPerPathway{posPathway}, i);
        end
    end
end

info = {};
for i = 1:length(uniquePathwayAll)
    info = [info;[uniquePathwayAll(i) {''} {''} {''} {''}]];
    for j = 1:length(posPerPathway{i});
        eqs1 = getRxn_cobraFormat(model, posPerPathway{i}(j),0);
        eqs2 = getRxn_cobraFormat(model, posPerPathway{i}(j),1);
        infoRxns = [model.rxns(posPerPathway{i}(j)) model.rxnNames(posPerPathway{i}(j)) eqs1 eqs2 model.grRules(posPerPathway{i}(j))];
        info = [info; infoRxns];
    end
    if i ~=length(uniquePathwayAll)
        info = [info;[{''} {''} {''} {''} {''}]];
    end
end

xlswrite(excelFile, info);

end