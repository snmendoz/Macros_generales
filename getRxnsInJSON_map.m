function rxns = getRxnsInJSON_map(file, geneSuffix)

fid = fopen(file, 'r+');
tline = fgetl(fid);
all = regexp(tline,'"bigg_id":"([a-zA-Z_0-9]+)"','tokens');
for i = 1:length(all)
    all{i} = all{i}{1};
end

genes = regexp(tline,['"bigg_id":"(' geneSuffix '[a-zA-Z_0-9]+)"'],'tokens');
for i = 1:length(genes)
    genes{i} = genes{i}{1};
end

mets = regexp(tline,'"bigg_id":"([a-zA-Z_0-9]+_[cpe])"','tokens');
for i = 1:length(mets)
    mets{i} = mets{i}{1};
end

rxns = setdiff(all, union(genes,mets));
fclose(fid);

end