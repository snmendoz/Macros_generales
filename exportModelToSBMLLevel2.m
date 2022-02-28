function exportModelToSBMLLevel2(model, modelName)

cd('D:\Dropbox\Databases\AuremeInputs');
fid = fopen('currentModels.txt','w+');
fprintf(fid,[modelName '\n']);
fclose(fid);

model = removeRxnsWithOutGeneAssociations(model);
% model.mets = regexprep(model.mets,'_((?!_).)+$','[$1]');
model.csense = repmat('E',length(model.mets),1);
file = [modelName '_rxnsOnlyGA'];
cd('D:\Dropbox\Databases\AuremeInputs')
if exist(modelName,'dir')~=7
    mkdir(modelName)
end
cd(modelName)

writeCbModel(model, 'fileName', file, 'format', 'sbml')
removeSquareBracketsFromBIGGModel([file '.xml'], 'metabolic_model_l3.sbml');

system('python "D:\\Dropbox\\Macros generales\\convertSBMLfroml3_to_l2.py"')
cd(['D:\Dropbox\Databases\AuremeInputs' filesep modelName]);
fixGeneAssociationInSBML('metabolic_model_cobra.xml', 'metabolic_model.sbml')
delete([modelName '_rxnsOnlyGA.xml']);
delete('metabolic_model_l3.sbml');
delete('metabolic_model_cobra.xml');

end