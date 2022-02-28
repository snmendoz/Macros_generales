load('D:\Dropbox\Databases\BIGG\bigg_84.mat')
bigg.csense = repmat('E', length(bigg.mets),1);
bigg = rmfield(bigg,'subSystems');
writeCbModel(bigg,'fileName', 'bigg_84.xml', 'format', 'sbml')