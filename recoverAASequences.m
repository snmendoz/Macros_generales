function sequences = recoverAASequences(model, folder, faaFileList, id)

sequences = cell(size(model.genes));
cd(folder);
is = 0;
for i = 1:length(faaFileList)
    fid = fopen(faaFileList{i});
    tline = fgetl(fid);
    while ischar(tline)
        
        if ~isempty(strfind(tline, '>'))
            locus_tag = regexp(tline, ['\[' id '=([^\]]*)\]'],'tokens');
            
            if ~isempty(locus_tag)
                locus_tag = regexprep(locus_tag{1}, '.*GeneID:','');
                [is, pos] = ismember(locus_tag, model.genes);
                if is
                    tline = fgetl(fid);
                    sequence = '';
                    while ischar(tline) && isempty(strfind(tline, '>'))
                        sequence = [sequence tline];
                        tline = fgetl(fid);
                    end
                    if isempty(sequences{pos})
                        sequences{pos} = sequence;
                    end
                end
            end
        end
        if is
            is = 0;
        else
            tline = fgetl(fid);
        end
        
    end
    fclose(fid);
end

end
