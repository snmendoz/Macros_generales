function exportCompoundstoAuReMe(model, pos, fileName)

fi = fopen(fileName, 'w+');

for i = 1:length(pos)
   metComp = regexp(model.mets{pos(i)},'.*\[(.*)\]$|.*_(.)$','tokens');
   metComp = metComp{1}{:};
   metName = regexprep(model.mets{pos(i)},'\[.*\]$|_.$','');
   fprintf(fi, '%s\t%s\n', metName, metComp);
end

fclose(fi);

end