function exportDataBaseToMeneco(database)

cd('D:\Dropbox\Databases\BIGG')
load(database)
exportBIGGModelToSBML(bigg,[database '.xml'],1)

ctsbml2_filepath = which('convertSBMLfroml3_to_l2_general.py');
baseFolder = 'D:\Dropbox\Databases\BIGG'; 
inputFileName = [database '.xml'];
outputFileName = [database '_l2.sbml'];
system(['python "' ctsbml2_filepath '" "' baseFolder '" "' inputFileName '" "' outputFileName '"'])

end