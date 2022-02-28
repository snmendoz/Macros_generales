function createTemplatesAuReMe(list)

cd('D:\Dropbox\Databases\AuremeInputs');
fid = fopen('currentModels.txt','w+');
for i = 1:length(list)
    fprintf(fid,[list{i} '\n']);
end
fclose(fid);

for i = 1:length(list)
    cd('D:\Dropbox\Databases\BIGG');
    load(list{i})
    model = eval(list{i});
    model = fixComparmentsInBIGGModel(model);
    model = removeRxnsWithOutGeneAssociations(model);
    model.mets = regexprep(model.mets,'_((?!_).)+$','[$1]');
    model.csense = repmat('E',length(model.mets),1);
    file = [list{i} '_rxnsOnlyGA'];
    cd('D:\Dropbox\Databases\AuremeInputs')
    if exist(list{i},'dir')~=7
        mkdir(list{i})
    end
    cd(list{i})
    if exist(['D:\Dropbox\Databases\BIGG\genomes\all\' list{i} '.faa'],'file') ==2
        copyfile(['D:\Dropbox\Databases\BIGG\genomes\all\' list{i} '.faa'], ...
           ['D:\Dropbox\Databases\AuremeInputs\' list{i} '\FAA_model.faa']);
    end
    
%     exportBIGGModelToSBML(model, file)
%     removeSquareBracketsFromBIGGModel([file '.xml'], 'metabolic_model_l3.sbml');
    
end

system('python "D:\\Dropbox\\Macros generales\\convertSBMLfroml3_to_l2.py"')

for i = 1:length(list)
    cd(['D:\Dropbox\Databases\AuremeInputs' filesep list{i}]);
    fixGeneAssociationInSBML('metabolic_model_cobra.xml', 'metabolic_model.sbml')
    delete([list{i} '_rxnsOnlyGA.xml']);
    delete('metabolic_model_l3.sbml');
    delete('metabolic_model_cobra.xml');
end

end