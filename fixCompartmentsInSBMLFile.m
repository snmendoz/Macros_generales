function fixCompartmentsInSBMLFile(fileName)

fid2 = fopen([fileName '_fixed.xml'],'w+');

fid = fopen([fileName '.xml'],'r+');
tline = fgetl(fid); 

while ischar(tline)
    if ~isempty(strfind(tline,'      <species metaid="'))
        metabolite = regexp(tline, 'species metaid="([^"]*)"','tokens');
        metabolite = metabolite{1};
        compartment = getCompartmentsFromMetList(metabolite);
        compartmentInSBML = regexp(tline, 'compartment="([^"]*)"','tokens');
        compartmentInSBML = compartmentInSBML{1};
        if strcmp(compartmentInSBML{1},'cell')
            
        end
        newLine = regexprep(tline,'compartment="cell"',['compartment="' compartment{1} '"']);
        fprintf(fid2, '%s\n', newLine);
    else
        fprintf(fid2, '%s\n', tline);
    end
    
    tline = fgetl(fid); 

end
fclose(fid);
fclose(fid2);


end


