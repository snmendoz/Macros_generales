function exportTargetsAuReMe(model, pos, fileName)

fi = fopen(fileName, 'w+');

for i = 1:length(pos)
   metName = model.mets{pos};
   comp = regexp(model.mets{pos},'.*\[(.*)\]$','tokens');
   fprintf(fi, '%s\t%s\n', metName, metComp);
end

fclose(fi);

end