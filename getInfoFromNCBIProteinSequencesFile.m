function info = getInfoFromNCBIProteinSequencesFile(fileName)

fi = fopen(fileName);
tline = fgetl(fi);

has_locus_tag = ~isempty(strfind(tline,'[locus_tag='));
has_protein = ~isempty(strfind(tline,'[protein='));
has_protein_id = ~isempty(strfind(tline,'[protein_id='));
has_location = ~isempty(strfind(tline,'[location='));
has_gbkey = ~isempty(strfind(tline,'[gbkey='));
% ids = cell(7000,1);
% sequences = cell(7000,1);
info.ids = cell(7000,1);
info.sequences = cell(7000,1);
if has_locus_tag; info.locus_tag = cell(7000,1); end;
if has_protein; info.protein = cell(7000,1); end;
if has_protein_id; info.protein_id = cell(7000,1); end;
if has_location; info.location = cell(7000,1); end;
if has_gbkey; info.gbkey = cell(7000,1); end;

cont = 0;
while ischar(tline)
    if ~isempty(strfind(tline, '>'))
        cont = cont+1;
        disp(cont)
        if has_locus_tag
            locus_tag = regexp(tline, '\[locus_tag=([^\[]*)\]','tokens');
            if ~isempty(locus_tag)
                info.locus_tag{cont}= locus_tag{1}{1};
            end
        end
        if has_protein
            protein = regexp(tline, '\[protein=([^\[]*)\]','tokens');
            if ~isempty(protein)
                info.protein{cont}= protein{1}{1};
            end
        end
        if has_protein_id
            protein_id = regexp(tline, '\[protein_id=([^\[]*)\]','tokens');
            if ~isempty(protein_id)
                info.protein_id{cont}= protein_id{1}{1};
            end
        end
        if has_location
            location = regexp(tline, '\[location=([^\[]*)\]','tokens');
            if ~isempty(location)
                info.location{cont}= location{1}{1};
            end
        end
        if has_gbkey
            gbkey = regexp(tline, '\[gbkey=([^\[]*)\]','tokens');
            if ~isempty(gbkey)
                info.gbkey{cont}= gbkey{1}{1};
            end
        end
        id = regexp(tline, '^>([^\[]*) [','tokens');
        info.ids(cont) =  id{1};
        
        tline = fgetl(fi);
        sequence = '';
        while ischar(tline) && isempty(strfind(tline, '>'))
            sequence = [sequence tline];
            tline = fgetl(fi);
        end
        info.sequences{cont} = sequence;
    end
end

info.ids = info.ids(1:cont);
info.sequences = info.sequences(1:cont);
if has_locus_tag; info.locus_tag = info.locus_tag(1:cont); end;
if has_protein; info.protein = info.protein(1:cont); end;
if has_protein_id; info.protein_id = info.protein_id(1:cont); end;
if has_location; info.location = info.location(1:cont); end;
if has_gbkey; info.gbkey = info.gbkey(1:cont); end;

fclose(fi);

end