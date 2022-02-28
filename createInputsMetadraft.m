function createInputsMetadraft

cd('D:\Dropbox\Databases\BIGG\genomes\all')
files = dir('*.gb');
for i = 1:length(files)
    copyfile(files(i).name,['D:\Dropbox\Databases\MetadraftInputs\' files(i).name]);
end

cd('D:\Dropbox\Databases\BIGG')
fi = fopen('modelNames.txt','r+');
for i = 1:length(files)
    name= fgetl(fi);
    copyfile([name '.xml'],['D:\Dropbox\Databases\MetadraftInputs\' name '.xml']);
end
fclose(fi);
end