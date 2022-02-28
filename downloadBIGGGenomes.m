function downloadBIGGGenomes
%open matlab using cmd
cd('D:\Dropbox\Macros generales')
system('python "D:\\Dropbox\\Macros generales\\getBIGGGeneBankFiles_general.py"')

s = tdfread('D:\Dropbox\Databases\BIGG\genomes\all\genomeIDs.txt');
modelNames = cellstr(s.ModelID);
genomeIDs = cellstr(s.GenomeID);
unGenomeIDs = unique(genomeIDs);

empty = find(strcmp(genomeIDs,'null'));
withOutGenome = modelNames(strcmp(genomeIDs,'null'));

cd('D:\Dropbox\Databases\BIGG\genomes\all')

for i = 1:length(modelNames)
    if exist([modelNames{i} '.faa'],'file')==2
        copyfile([modelNames{i} '.faa'], [modelNames{i} '_copy.faa']);
        delete([modelNames{i} '.faa']);
    end
    if exist([modelNames{i} '.gb'],'file')==2
        TransformarFASTA_ParaPantographGeneral([modelNames{i} '.gb'],[modelNames{i} '.faa']);
    end
end

end
