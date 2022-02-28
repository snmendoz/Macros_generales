function renameMetadraftModels

cd('D:\Dropbox\new_templates\bigg2')

files = dir('*.seqplus.xml');
for i = 1:length(files)
    file_i = files(i).name;
    if isempty(strfind(file_i,'usr-'))
        continue;
    end
    
    pos1 = strfind(file_i, '-'); pos1 = pos1(1);
    pos2 = strfind(file_i, ')-');
    modelName = file_i(pos1+1:pos2-1);
    
    f = fopen([modelName '.gb'],'r+');
    
    tline = fgetl(f);
    while ischar(tline)
        if ~isempty(strfind(tline, 'ORGANISM'))
            pos1 = strfind(tline, 'ORGANISM');
            speciesName = tline(pos1+10:end);
            posEmpty = strfind(speciesName, ' ');
            posStr = strfind(speciesName, ' str.');
            posSubStr = strfind(speciesName, ' substr.');
            posSubSpe = strfind(speciesName, ' subsp.');
            abbr = [speciesName(1) speciesName(posEmpty+1:posEmpty+2)];
            %             if ~isempty(posSubSpe)
            %                 first = posEmpty(find(posEmpty>posSubSpe));
            %                 if length(first)>1
            %                     first = first(2);
            %                     abbr = [abbr '_' speciesName(posSubSpe+8:first-1)];
            %                 else
            %                     abbr = [abbr '_' speciesName(posSubSpe+8:end)];
            %                 end
            %             end
            %
            %             if ~isempty(posStr)
            %                 first = posEmpty(find(posEmpty>posStr));
            %                 if length(first)>1
            %                     first = first(2);
            %                     abbr = [abbr '_' speciesName(posStr+6:first-1)];
            %                 else
            %                     abbr = [abbr '_' speciesName(posStr+6:end)];
            %                 end
            %             end
            %
            %             if ~isempty(posSubStr)
            %                 first = posEmpty(find(posEmpty>posSubStr));
            %                 if length(first)>1
            %                     first = first(2);
            %                     abbr = [abbr '_' speciesName(posSubStr+9:first-1)];
            %                 else
            %                     abbr = [abbr '_' speciesName(posSubStr+9:end)];
            %                 end
            %             end
            
            abbr = regexprep(lower(abbr),{'-','/'},{'',''});
            %             disp(speciesName)
            disp(abbr);
            %             disp(modelName)
            break;
        end
        tline = fgetl(f);
    end
    fclose(f);
    value = abbr;
    
    copyfile(file_i, regexprep(file_i, ['usr-',modelName],['bigg2-',value]));
    delete(file_i);
    
end

files = dir('*.seqplus.json');
for i = 1:length(files)
    file_i = files(i).name;
    if isempty(strfind(file_i,'usr-'))
        continue;
    end
    pos1 = strfind(file_i, '-'); pos1 = pos1(1);
    pos2 = strfind(file_i, ')-');
    modelName = file_i(pos1+1:pos2-1);
    
    f = fopen([modelName '.gb'],'r+');
    
    tline = fgetl(f);
    while ischar(tline)
        if ~isempty(strfind(tline, 'ORGANISM'))
            pos1 = strfind(tline, 'ORGANISM');
            speciesName = tline(pos1+10:end);
            posEmpty = strfind(speciesName, ' ');
            posStr = strfind(speciesName, ' str.');
            posSubStr = strfind(speciesName, ' substr.');
            posSubSpe = strfind(speciesName, ' subsp.');
            abbr = [speciesName(1) speciesName(posEmpty+1:posEmpty+2)];
            %             if ~isempty(posSubSpe)
            %                 first = posEmpty(find(posEmpty>posSubSpe));
            %                 if length(first)>1
            %                     first = first(2);
            %                     abbr = [abbr '_' speciesName(posSubSpe+8:first-1)];
            %                 else
            %                     abbr = [abbr '_' speciesName(posSubSpe+8:end)];
            %                 end
            %             end
            %
            %             if ~isempty(posStr)
            %                 first = posEmpty(find(posEmpty>posStr));
            %                 if length(first)>1
            %                     first = first(2);
            %                     abbr = [abbr '_' speciesName(posStr+6:first-1)];
            %                 else
            %                     abbr = [abbr '_' speciesName(posStr+6:end)];
            %                 end
            %             end
            %
            %             if ~isempty(posSubStr)
            %                 first = posEmpty(find(posEmpty>posSubStr));
            %                 if length(first)>1
            %                     first = first(2);
            %                     abbr = [abbr '_' speciesName(posSubStr+9:first-1)];
            %                 else
            %                     abbr = [abbr '_' speciesName(posSubStr+9:end)];
            %                 end
            %             end
            
            abbr = regexprep(lower(abbr),{'-','/'},{'',''});
            %             disp(speciesName)
            disp(abbr);
            %             disp(modelName)
            break;
        end
        tline = fgetl(f);
    end
    fclose(f);
    value = abbr;
    
    copyfile(file_i, regexprep(file_i, ['usr-',modelName],['bigg2-',value]));
    delete(file_i);
    
end

end