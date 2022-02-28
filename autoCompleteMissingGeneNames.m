function model = autoCompleteMissingGeneNames(model)

for i = 1:length(model.geneNames)
    if isempty(model.geneNames{i})
        model.geneNames{i} = model.genes{i};
    end
end

end