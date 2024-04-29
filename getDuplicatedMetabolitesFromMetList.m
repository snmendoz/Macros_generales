function [duplicatedMets, posMets, posAll] = getDuplicatedMetabolitesFromMetList(metList)

[~, pos ] = unique(metList); 
posMets = setdiff(1:length(metList), pos);
duplicatedMets = metList(posMets);
posAll = cellfun(@(x) find(strcmp(metList, x)), duplicatedMets, 'UniformOutput', false);
end