function exportBoundsToCSV(model, fileName)

info = [model.rxns, num2cell(model.lb) num2cell(model.ub)];
for i = 1:size(info,1)
    for j = 2:size(info,2)
        info{i,j} = num2str(info{i,j});
    end
end
exportToCSV(fileName, info)

end