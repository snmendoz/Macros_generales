function WriteReactionsToRemoveInAuremeFormat(fileName, forAdding)

%agregadas desde bigg
fi=fopen([ fileName '.txt'],'w+');
fprintf(fi,'idRef\tComment\tAction\n');
for i = 1:length(forAdding)
    fprintf(fi,[forAdding{i} '\tmanual curation\tdelete\n']);
end
fclose(fi);

end