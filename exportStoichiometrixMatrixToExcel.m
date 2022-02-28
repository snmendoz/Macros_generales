function exportStoichiometrixMatrixToExcel(model, fileName)

matrix = full(model.S);
info = [{''} model.rxns'; [model.mets, num2cell(matrix)]];
xlswrite(fileName,info, 'S')

end