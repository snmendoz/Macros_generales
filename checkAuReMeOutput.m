function [modelIsOK, rxnDiff, metDiff, same] = checkAuReMeOutput(folder, file)

load('D:\Dropbox\Databases\BIGG\bigg_85.mat')

cd(folder)
model = readCbModel(file);
model.mets = regexprep(model.mets,'\[(.)\]','_$1');

rxnDiff = setdiff(model.rxns,bigg.rxns);
metDiff = setdiff(model.mets,bigg.mets);

same = zeros(size(model.rxns));
for i = 1:length(model.rxns)
    [is, pos]= ismember(model.rxns(i), bigg.rxns);
    if is
        [sameRxnID, sameRxnFormula, writtenInSameSense, sameReversibility] = compareRxns(model, i, bigg, pos);
        if sameRxnFormula
            same(i) = 1;
        end
    end
end

modelIsOK = 0;
if isempty(rxnDiff) && isempty(metDiff) && isempty(find(same==0))
    modelIsOK = 1;
end

end