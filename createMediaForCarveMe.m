function createMediaForCarveMe(mediaID, mediaName, compoundIDs,  compoundNames)

fileID = fopen([mediaID '.tsv'], 'w+');

for i = 1:length(compoundIDs)
    fprintf(fileID, '%s\t%s\t%s\t%s\n', mediaID, mediaName, compoundIDs{i}, compoundNames{i});
end

fclose(fileID);

end