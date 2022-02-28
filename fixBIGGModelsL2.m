function fixBIGGModelsL2

cd('D:\Dropbox\Databases\BIGG');
fid = fopen('modelNames.txt');

cd('D:\Dropbox\Databases\BIGG\SBML_L2');
tline = fgetl(fid);
while ischar(tline)
    if exist([tline '_l2.xml'],'file')==2 && exist([tline '_l2_v1.xml'],'file')~=2
        fixGeneAssociationInSBML([tline '_l2.xml'], [tline '_l2_v1.xml']);
    end
    tline = fgetl(fid);
end

fclose(fid);

end

