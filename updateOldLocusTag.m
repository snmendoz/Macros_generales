function model = updateOldLocusTag(model, infoGenes)

if isfield(model, 'geneOldLocusTag')
    for i = 1:length(model.genes)
        pos = find(strcmp(model.genes{i},infoGenes.locus_tag));
        if ~isempty(infoGenes.old_locus_tag{pos})
            model.geneOldLocusTag{i} = infoGenes.old_locus_tag{pos};
        end
    end
end
end