function fixGeneAssociationInSBML(inputFile, outputFile)

fID = fopen(inputFile, 'r+');
fID2 = fopen(outputFile, 'w+');
tline = fgetl(fID);
while ischar(tline)
    if ~isempty(strfind(tline, '<p>GENE_ASSOCIATION'))
        str = tline;
        str = regexprep(str, 'GENE_ASSOCIATION :', 'GENE_ASSOCIATION:');
        fprintf(fID2, [str '\n']);
    else
        fprintf(fID2, [tline '\n']);
        
    end
    tline = fgetl(fID);
end

fclose(fID);
fclose(fID2);

end