function exportSolutionToJson(model,solution, fileName)

if ~endsWith(fileName, '.json')
    fileName = [fileName '.json'];
end

s = arrayfun(@(x) num2str(x), solution, 'UniformOutput', false);

lstring = regexprep(strjoin(strcat('"', model.rxns, '":', arrayfun(@(x) num2str(x), solution, 'UniformOutput', false)),', '), ':',': ');

string = ['{' lstring '}'];

fid = fopen(fileName,'w+');
fprintf(fid, string);

fclose(fid);

end