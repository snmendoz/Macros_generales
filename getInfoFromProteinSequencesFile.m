function info = getInfoFromProteinSequencesFile(fileName)

fi = fopen(fileName);
tline = fgetl(fi);

info.ids = cell(7000,1);
info.sequences = cell(7000,1);
cont = 0;
while ischar(tline)
    if ~isempty(strfind(tline, '>'))
        cont = cont+1;
        disp(cont)
        locus_tag = tline(2:end);
        info.ids{cont}= locus_tag;
    end
    tline = fgetl(fi);
    sequence = '';
    while ischar(tline) && isempty(strfind(tline, '>'))
        sequence = [sequence tline];
        tline = fgetl(fi);
    end
    info.sequences{cont} = sequence;
end

info.ids = info.ids(1:cont);
info.sequences = info.sequences(1:cont);
fclose(fi);



end