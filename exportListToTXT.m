function exportListToTXT(list,fileName)

if ~isempty(strfind(fileName,'.txt'))
    fileID = fopen(fileName,'wt');
else
    fileID = fopen([fileName '.txt'],'wt');
end
for i=1:length(list)
    fprintf(fileID, '%s\n', list{i});
end
fclose(fileID);

end