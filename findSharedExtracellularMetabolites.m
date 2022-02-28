function shared = findSharedExtracellularMetabolites(models)

ExtMetModel1 = findExtracellularMetabolites(models{1});
ExtMetModel2 = findExtracellularMetabolites(models{2});

% ExtMetModel1 = intersect(findExtracellularMetabolites(models{1}), regexprep(models{1}.rxns(findExcRxnsWithIDs(models{1})), 'EX_',''));
% ExtMetModel2 = intersect(findExtracellularMetabolites(models{2}), regexprep(models{2}.rxns(findExcRxnsWithIDs(models{2})), 'EX_',''));

shared = intersect(ExtMetModel1, ExtMetModel2);

if length(models)>2
    for i = 3:length(models)
        ExtMetModeln = findExtracellularMetabolites(models{i});
        shared = intersect(shared, ExtMetModeln);
    end
end


end