function [rxnsMedia, valuesMedia] = getMediumFromExcelFile(fileName, sheet)

[~, data] = xlsread(fileName,sheet);
basalMedium_rxns = data(1:end,2);
rxnsMedia = {};
for k = 1:length(basalMedium_rxns)
    rxnsMedia = union(rxnsMedia, strsplit(basalMedium_rxns{k},';'));
end
rxnsMedia = setdiff(rxnsMedia,'');
valuesMedia = ones(size(rxnsMedia))*10;

end