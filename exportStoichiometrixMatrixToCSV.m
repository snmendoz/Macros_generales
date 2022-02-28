function exportStoichiometrixMatrixToCSV(model, fileName)

matrix = full(model.S);
info = [{''} model.rxns'; [model.mets, num2cell(matrix)]];
for i = 2:size(info,1)
    for j = 2:size(info,2)
        info{i,j} = num2str(info{i,j});
    end
end
exportToCSV(fileName, info)

end