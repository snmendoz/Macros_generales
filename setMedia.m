function [model, rxnsMedia, compounds, rxnsPerCompound]= setMedia(model, excelFile, sheet, maxUptakeRate)

if nargin < 4 || isempty(maxUptakeRate)
    maxUptakeRate = 10;
end

[~, data] = xlsread(excelFile,sheet);
basalMedium = data(1:end,1);
compounds = data(1:end,2);
rxnsMedia = {};
for i = 1:length(basalMedium)
    rxnsMedia = union(rxnsMedia, strsplit(compounds{i},';'));
    rxnsPerCompound{i} = strcat('EX_',setdiff(strsplit(compounds{i},';'),''),'_e');
end
rxnsMedia = strcat('EX_',setdiff(rxnsMedia,''),'_e');
lowerValuesBasalMedium = -ones(size(rxnsMedia))*maxUptakeRate;

pos_Ex = find(cellfun(@isempty, strfind(model.rxns,'EX_'))==0);
model = changeRxnBounds(model, model.rxns(pos_Ex), 0, 'l');
model = changeRxnBounds(model, model.rxns(pos_Ex(find(model.ub(pos_Ex)<0))), 0, 'u');

model = changeRxnBounds(model, rxnsMedia, lowerValuesBasalMedium, 'l');

end