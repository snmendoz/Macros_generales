function model = addMNXIdentifiersFromSourceDatabase(model,listOfMetFieldPriorities, listOfRxnFieldPriorities, overwrite, ...
    metMNXxrefSourceID, metMNXxrefSourceMet, metMNXIDs, rxnMNXxrefSourceID, rxnMNXxrefSourceRxn, rxnMNXxrefID)

metDatabases = {'biggM','keggC','keggD','keggE','keggG','metacycM','metacycM','seedM'};
rxnDatabases = {'biggR','keggR','keggR','keggR','keggR','metacycR','metacycR','seedR'};
fields = {'BiGGID','KEGGID','KEGGID','KEGGID','KEGGID','MetaCycID','BioCycID','SEEDID'};

metMetaNetXIDs = repmat({''},length(model.mets),1);

for i = 1:length(listOfMetFieldPriorities)
    fieldValues =  eval(['model.' listOfMetFieldPriorities{i}]);
    pos_notEmpty = intersect(find(cellfun(@isempty, fieldValues)==0), find(cellfun(@isempty, metMetaNetXIDs)==1));
    met_ids = fieldValues(pos_notEmpty);
    pos = getPosOfElementsInArray(listOfMetFieldPriorities(i),strcat('met',fields));
    for k = 1:length(pos)
        inputLanguage = metDatabases{pos(k)};
        metIDs_outputLanguage = getMetIDsInTargetLanguageFromInputLanguage_4_0(met_ids,...
                inputLanguage,'MNX', metMNXxrefSourceID, metMNXxrefSourceMet, metMNXIDs);
        for j = 1:length(met_ids)
            if isempty(metMetaNetXIDs{pos_notEmpty(j)}) && (isempty(model.metMetaNetXID{pos_notEmpty(j)}) || (~isempty(model.metMetaNetXID{pos_notEmpty(j)}) && overwrite) )
                metMetaNetXIDs{pos_notEmpty(j)} = metIDs_outputLanguage{j};
                model.metMetaNetXID{pos_notEmpty(j)} = metIDs_outputLanguage{j};
            end
        end
    end
end

model.metMetaNetXID = metMetaNetXIDs;


rxnMetaNetXIDs = repmat({''},length(model.rxns),1);

for i = 1:length(listOfRxnFieldPriorities)
    fieldValues =  eval(['model.' listOfRxnFieldPriorities{i}]);
    pos_notEmpty = intersect(find(cellfun(@isempty, fieldValues)==0), find(cellfun(@isempty, rxnMetaNetXIDs)==1));
    rxn_ids = fieldValues(pos_notEmpty);
    pos = getPosOfElementsInArray(listOfRxnFieldPriorities(i),strcat('rxn',fields));
    for k = 1:length(pos)
        inputLanguage = rxnDatabases{pos(k)};
        rxnIDs_outputLanguage = getRxnIDsInTargetLanguageFromInputLanguage_4_0(rxn_ids,...
                inputLanguage,'MNX', rxnMNXxrefSourceID, rxnMNXxrefSourceRxn, rxnMNXxrefID);
        for j = 1:length(rxn_ids)
            if isempty(rxnMetaNetXIDs{pos_notEmpty(j)}) && (isempty(model.rxnMetaNetXID{pos_notEmpty(j)}) || (~isempty(model.rxnMetaNetXID{pos_notEmpty(j)}) && overwrite) )
                rxnMetaNetXIDs{pos_notEmpty(j)} = rxnIDs_outputLanguage{j};
                model.rxnMetaNetXID{pos_notEmpty(j)} = rxnIDs_outputLanguage{j};
            end
        end
    end
end

model.rxnMetaNetXID = rxnMetaNetXIDs;

end