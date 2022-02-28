function media = readMediaFromTXT(fileName, separator)

info = readCSVFile_general(fileName, separator);
rxns = info(2:end,1);
lowerBounds = cellfun(@ str2num,info(2:end,2),'UniformOutput',0);

media.reactions = rxns;
media.lb = cell2mat(lowerBounds);

end