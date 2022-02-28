function exportOnlySequencesToFastaFile(info, fileName)

fid = fopen(fileName,'w+');
for i = 1:length(info.ids)
    fprintf(fid, [info.sequences{i} '\n']);
end
fclose(fid);

end