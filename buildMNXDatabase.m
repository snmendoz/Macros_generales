function buildMNXDatabase(version)

cd(['D:\Dropbox\Databases\MNX\' version])


% %% metabolite info
% fid = fopen('chem_prop.tsv');
% tline = fgetl(fid);
% while ischar(tline)
%     if startsWith(tline,'#ID')
%         break;
%     end
%     tline = fgetl(fid);
% end
% tline = fgetl(fid);
% tline = fgetl(fid);
% fields = split(regexprep(tline,'\t','^'),'^');
% cont = 0;
% info = repmat({''},1100000,length(fields));
% 
% while ischar(tline) 
%     disp(cont)
%     fields = strsplit(tline,'\t');
%     if length(fields)<length(fields)
%         fields = split(regexprep(tline,'\t','^'),'^');
%     end
%     if length(fields)~=length(fields)
%         fields = split(regexprep(tline,'\t','^'),'^');
%     end
%     cont = cont+1;
%     info(cont,:) = fields';
%     tline = fgetl(fid);
% end
% 
% info = info(1:cont,:);
% metMNXID = info(:,1);
% metMNXname = info(:,2);
% metMNXreference = info(:,3);
% metMNXformula = info(:,4);
% metMNXcharge = info(:,5);
% metMNXmass = info(:,6);
% metMNXInChI = info(:,7);
% metMNXInChIKey = info(:,8);
% metMNXSMILE = info(:,9);
% 
% save('metMNXprop_4_0','info')
% save('metMNXID_4_0','metMNXID')
% save('metMNXname_4_0','metMNXname')
% save('metMNXreference_4_0','metMNXreference')
% save('metMNXformula_4_0','metMNXformula')
% save('metMNXcharge_4_0','metMNXcharge')
% save('metMNXmass_4_0','metMNXmass')
% save('metMNXInChI_4_0','metMNXInChI')
% save('metMNXInChIKey_4_0','metMNXInChIKey')
% save('metMNXSMILE_4_0','metMNXSMILE')
% 
% fclose(fid);
%% metabolite mapping
fid = fopen('chem_xref.tsv');
tline = fgetl(fid);
while ischar(tline)
    if startsWith(tline,'#source')
        break;
    end
    tline = fgetl(fid);
end
tline = fgetl(fid);
tline = fgetl(fid);
fields = split(regexprep(tline,'\t','^'),'^');
cont = 0;
info = repmat({''},1400000,length(fields));

while ischar(tline) 
    disp(cont)
    fields = strsplit(tline,'\t');
    if length(fields)<length(fields)
        fields = split(regexprep(tline,'\t','^'),'^');
    end
    if length(fields)~=length(fields)
        fields = split(regexprep(tline,'\t','^'),'^');
    end
    cont = cont+1;
    info(cont,:) = fields';
    tline = fgetl(fid);
end

info = info(1:cont,:);
metMNXxrefSource = info(:,1);
metMNXxrefID = info(:,2);

posToRemove = find(cellfun(@(x) ismember(':', x), metMNXxrefSource)==0);
metMNXxrefID(posToRemove) = [];
metMNXxrefSource(posToRemove) = [];

metMNXxrefSourceID = repmat({''},length(metMNXxrefSource),1);
metMNXxrefSourceMet = repmat({''},length(metMNXxrefSource),1);

for i = 1:length(metMNXxrefSource)
    split_data_i = regexp(metMNXxrefSource{i}, '([^:]*):(.*)','tokens');
    metMNXxrefSourceID{i} = split_data_i{1}{1};
    metMNXxrefSourceMet{i} = split_data_i{1}{2};
end

% split_data = split(metMNXxrefSource,':');
% metMNXxrefSourceID = split_data(:,1);
% metMNXxrefSourcerxn = split_data(:,2);

save('metMNXxref_4_0','info')
save('metMNXxrefSource_4_0','metMNXxrefSource')
save('metMNXxrefID_4_0','metMNXxrefID')
save('metMNXxrefSourceID_4_0','metMNXxrefSourceID')
save('metMNXxrefSourceMet_4_0','metMNXxrefSourceMet')


fclose(fid);

%% reactions
fid2 = fopen('reac_prop.tsv');
tline = fgetl(fid2);
while ischar(tline)
    if startsWith(tline,'#ID')
        break;
    end
    tline = fgetl(fid2);
end
tline = fgetl(fid2);
tline = fgetl(fid2);
fields = split(regexprep(tline,'\t','^'),'^');
cont = 0;
info = repmat({''},100000,length(fields));

while ischar(tline) 
    disp(cont)
    fields = strsplit(tline,'\t');
    if length(fields)<6
        fields = split(regexprep(tline,'\t','^'),'^');
    end
    if length(fields)~=6
        fields = split(regexprep(tline,'\t','^'),'^');
    end
    cont = cont+1;
    info(cont,:) = fields';
    tline = fgetl(fid2);
end

info = info(1:cont,:);
rxnMNXID = info(:,1);
rxnMNXequation = info(:,2);
rxnMNXreference = info(:,3);
rxnMNXclassifs = info(:,4);
rxnMNXis_balanced = info(:,5);
rxnMNXis_transport = info(:,6);

save('rxnMNXprop','info')
save('rxnMNXID','rxnMNXID')
save('rxnMNXequation','rxnMNXequation')
save('rxnMNXreference','rxnMNXreference')
save('rxnMNXclassifs','rxnMNXclassifs')
save('rxnMNXis_balanced','rxnMNXis_balanced')
save('rxnMNXis_transport','rxnMNXis_transport')

fclose(fid2);


%% reaction mapping
fid = fopen('reac_xref.tsv');
tline = fgetl(fid);
while ischar(tline)
    if startsWith(tline,'#source')
        break;
    end
    tline = fgetl(fid);
end
tline = fgetl(fid);
tline = fgetl(fid);
fields = split(regexprep(tline,'\t','^'),'^');
cont = 0;
info = repmat({''},200000,length(fields));

while ischar(tline) 
    disp(cont)
    fields = strsplit(tline,'\t');
    if length(fields)<length(fields)
        fields = split(regexprep(tline,'\t','^'),'^');
    end
    if length(fields)~=length(fields)
        fields = split(regexprep(tline,'\t','^'),'^');
    end
    cont = cont+1;
    info(cont,:) = fields';
    tline = fgetl(fid);
end

info = info(1:cont,:);
rxnMNXxrefSource = info(:,1);
rxnMNXxrefID = info(:,2);

posToRemove = find(cellfun(@(x) ismember(':', x), rxnMNXxrefSource)==0);
rxnMNXxrefID(posToRemove) = [];
rxnMNXxrefSource(posToRemove) = [];

rxnMNXxrefSourceID = repmat({''},length(rxnMNXxrefSource),1);
rxnMNXxrefSourceRxn = repmat({''},length(rxnMNXxrefSource),1);

for i = 1:length(rxnMNXxrefSource)
    split_data_i = regexp(rxnMNXxrefSource{i}, '([^:]*):(.*)','tokens');
    rxnMNXxrefSourceID{i} = split_data_i{1}{1};
    rxnMNXxrefSourceRxn{i} = split_data_i{1}{2};
end


% split_data = split(rxnMNXxrefSource,':');
% rxnMNXxrefSourceID = split_data(:,1);
% rxnMNXxrefSourceRxn = split_data(:,2);

save('rxnMNXxref_4_0','info')
save('rxnMNXxrefSource_4_0','rxnMNXxrefSource')
save('rxnMNXxrefID_4_0','rxnMNXxrefID')
save('rxnMNXxrefSourceID_4_0','rxnMNXxrefSourceID')
save('rxnMNXxrefSourceRxn_4_0','rxnMNXxrefSourceRxn')
fclose(fid);

end