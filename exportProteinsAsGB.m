function exportProteinsAsGB(info, fileName)

fid = fopen(fileName,'w+');
fprintf(fid, 'FEATURES             Location/Qualifiers\n');

for i = 1:length(info.genes)
    fprintf(fid, ['     gene            ' info.locations{i} '\n']);
    fprintf(fid, ['                     /locus_tag="' info.genes{i} '"\n']);
    fprintf(fid, ['                     /db_xref="' info.geneIDs{i} '"\n']);
    fprintf(fid, ['     CDS             ' info.locations{i} '\n']);
    fprintf(fid, ['                     /locus_tag="' info.genes{i} '"\n']);
    fprintf(fid, ['                     /db_xref="' info.geneIDs{i} '"\n']);
    if length(info.sequences{i})<=44
        fprintf(fid, ['                     /translation="' info.sequences{i} '"\n']);
    else
        fprintf(fid, ['                     /translation="' info.sequences{i}(1:44) '\n']);
        rest = info.sequences{i}(45:end);
        a = floor(length(rest)/58);
        for j = 1:a
            fprintf(fid, ['                     ' rest(58*(j-1)+1 : (58*(j-1)+1)+57) '\n']);
        end
        j = j+1;
        fprintf(fid, ['                     ' rest(58*(j-1)+1 : end) '"\n']);
    end
end

fprintf(fid, 'ORIGIN      ');
fclose(fid);


end