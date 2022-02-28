function exportMediaForCarveMe(mediaID, mediaName, mediaRxns, mediaFile)

fid = fopen(mediaFile,'w+');

mediaCompounds = regexprep(mediaRxns, {'^EX_','_e$'}, {'',''});

for i = 1:length(mediaRxns)
    fprintf(fid, '%s\t%s\t%s\t%s\n',mediaID, mediaName, mediaCompounds{i},mediaCompounds{i});
end

end