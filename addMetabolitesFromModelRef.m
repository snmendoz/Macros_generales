function modelOut = addMetabolitesFromModelRef(model, metaboliteIDs, modelRef)

modelOut = model;
overwrite = 1;
modelOut = addMetabolite(modelOut,metaboliteIDs);
posMets = getPosOfElementsInArray(metaboliteIDs, modelOut.mets);

metWOCompsReference = removeCompartmentFromMets(modelRef.mets);
metWOCompsUserModel = removeCompartmentFromMets(modelOut.mets);

modelFields = fieldnames(modelOut);
metFields_cells = modelFields(intersect(find(cellfun(@isempty, regexp(modelFields,'^met'))==0),find(cellfun(@(x) iscell(modelOut.(x)),modelFields))));
metFields_cells = setdiff(metFields_cells,'mets');
metFields_double = modelFields(intersect(find(cellfun(@isempty, regexp(modelFields,'^met'))==0),find(cellfun(@(x) isnumeric(modelOut.(x)),modelFields))));


if allMetsInCBMPYFormat(modelOut.mets) && allMetsInCOBRAFormat(modelRef.mets)
    modelRef = transformModelToCBMPYFormat(modelRef);
end
if allMetsInCBMPYFormat(modelRef.mets) && allMetsInCOBRAFormat(modelOut.mets)
    modelRef = transformModelToCOBRAFormat(modelRef);
end

for i = 1:length(posMets);
    posMetInModel_i = posMets(i);
    if ismember(metWOCompsUserModel(posMetInModel_i), metWOCompsReference)
        posMetInReference_i = find(strcmp(metWOCompsReference,metWOCompsUserModel{posMetInModel_i}));
        if length(posMetInReference_i)>1;
            posMetInReference_i = posMetInReference_i(1);
        end
        for j = 1:length(metFields_cells)
            if isfield(modelRef,metFields_cells{j}) && isfield(modelOut,metFields_cells{j}) && ~isempty(modelRef.(metFields_cells{j}){posMetInReference_i}) && (isempty(modelOut.(metFields_cells{j}){posMetInModel_i}) || (~isempty(modelOut.(metFields_cells{j}){posMetInModel_i}) && overwrite ))
                modelOut.(metFields_cells{j}){posMetInModel_i} = regexprep(modelRef.(metFields_cells{j}){posMetInReference_i},',',';');
            end
        end
        
        for j = 1:length(metFields_double)
            if isfield(modelRef,metFields_double{j}) && isfield(modelOut,metFields_double{j}) && ~isempty(modelRef.(metFields_double{j})(posMetInReference_i)) && (isempty(modelOut.(metFields_double{j})(posMetInModel_i)) || (~isempty(modelOut.(metFields_double{j})(posMetInModel_i)) && overwrite ))
                modelOut.(metFields_double{j})(posMetInModel_i) = modelRef.(metFields_double{j})(posMetInReference_i);
            end
        end
    end
end


end