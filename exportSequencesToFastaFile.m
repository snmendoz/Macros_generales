function exportSequencesToFastaFile(info, fileName)

fid = fopen(fileName,'w+');
for i = 1:length(info.ids)
    fprintf(fid, ['>' info.ids{i} '\n']);
    fprintf(fid, [info.sequences{i} '\n']);
end
fclose(fid);

end