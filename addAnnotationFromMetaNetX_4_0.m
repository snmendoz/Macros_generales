function modelOut = addAnnotationFromMetaNetX_4_0(model,  posRxns, posMets, overwrite, rxnMNXxrefID, rxnMNXxrefSourceID, rxnMNXxrefSourceRxn, metMNXxrefID, metMNXxrefSourceID, metMNXxrefSourceMet)

global rootFolder

if nargin<2 || isempty(posRxns); posRxns = 1:length(model.rxns); end
if nargin<3 || isempty(posMets); posMets = 1:length(model.mets); end; 
if nargin<4 || isempty(overwrite); overwrite = 0; end;
if nargin<5
    load(fullfile(rootFolder, 'MNX\rxnMNXxrefID'));
    load(fullfile(rootFolder, 'MNX\rxnOtherIDs'));
    load(fullfile(rootFolder, 'MNX\metMNXIDs.mat'));
    load(fullfile(rootFolder, 'MNX\metOtherIDs.mat'));
end

modelOut = model;

for i = 1:length(posMets)
    posMetInModel_i = posMets(i);
    if ~isempty(modelOut.metMetaNetXID{posMetInModel_i})
        id = modelOut.metMetaNetXID(posMetInModel_i);
        biggID = getMetIDsInTargetLanguageFromInputLanguage_4_0(id,'MNX','biggM', metMNXxrefSourceID, metMNXxrefSourceMet, metMNXxrefID);
        if ~isempty(biggID{1}) && ((~isempty(modelOut.metBiGGID{posMetInModel_i}) && overwrite) || isempty(modelOut.metBiGGID{posMetInModel_i})); modelOut.metBiGGID{posMetInModel_i} = strjoin(biggID{1}, ';'); end;
        keggID =  getMetIDsInTargetLanguageFromInputLanguage_4_0(id,'MNX','keggC', metMNXxrefSourceID, metMNXxrefSourceMet, metMNXxrefID);
        if ~isempty(keggID{1}) && ((~isempty(modelOut.metKEGGID{posMetInModel_i}) && overwrite) || isempty(modelOut.metKEGGID{posMetInModel_i})); modelOut.metKEGGID{posMetInModel_i} = strjoin(keggID{1}, ';'); end;
        seedID = getMetIDsInTargetLanguageFromInputLanguage_4_0(id,'MNX','seedM', metMNXxrefSourceID, metMNXxrefSourceMet, metMNXxrefID);
        if ~isempty(seedID{1}) && ((~isempty(modelOut.metSEEDID{posMetInModel_i}) && overwrite) || isempty(modelOut.metSEEDID{posMetInModel_i})); modelOut.metSEEDID{posMetInModel_i} = strjoin(seedID{1}, ';'); end;
        metacycID = getMetIDsInTargetLanguageFromInputLanguage_4_0(id,'MNX','metacycM', metMNXxrefSourceID, metMNXxrefSourceMet, metMNXxrefID);
        if ~isempty(metacycID{1}) && ((~isempty(modelOut.metMetaCycID{posMetInModel_i}) && overwrite) || isempty(modelOut.metMetaCycID{posMetInModel_i})); modelOut.metMetaCycID{posMetInModel_i} = strjoin(metacycID{1}, ';'); end;
        chebiID = getMetIDsInTargetLanguageFromInputLanguage(id,'MNX','chebi', metMNXxrefSourceID, metMNXxrefSourceMet, metMNXIDs);
        if ~isempty(chebiID{1}) && ((~isempty(modelOut.metChEBIID{posMetInModel_i}) && overwrite) || isempty(modelOut.metChEBIID{posMetInModel_i})); modelOut.metChEBIID{posMetInModel_i} = strjoin(chebiID{1}, ';'); end;
        envipathID = getMetIDsInTargetLanguageFromInputLanguage(id,'MNX','envipath', metMNXxrefSourceID, metMNXxrefSourceMet, metMNXIDs);
        if ~isempty(envipathID{1}) && ((~isempty(modelOut.metEnvipathID{posMetInModel_i}) && overwrite) || isempty(modelOut.metEnvipathID{posMetInModel_i})); modelOut.metEnvipathID{posMetInModel_i} = strjoin(envipathID{1}, ';'); end;
        hmdbID = getMetIDsInTargetLanguageFromInputLanguage(id,'MNX','hmdb', metMNXxrefSourceID, metMNXxrefSourceMet, metMNXIDs);
        if ~isempty(hmdbID{1}) && ((~isempty(modelOut.metHMDBID{posMetInModel_i}) && overwrite) || isempty(modelOut.metHMDBID{posMetInModel_i})); modelOut.metHMDBID{posMetInModel_i} = strjoin(hmdbID{1}, ';'); end;
        lipidmapsID = getMetIDsInTargetLanguageFromInputLanguage(id,'MNX','lipidmaps', metMNXxrefSourceID, metMNXxrefSourceMet, metMNXIDs);
        if ~isempty(lipidmapsID{1}) && ((~isempty(modelOut.metLipidMapsID{posMetInModel_i}) && overwrite) || isempty(modelOut.metLipidMapsID{posMetInModel_i})); modelOut.metLipidMapsID{posMetInModel_i} = strjoin(lipidmapsID{1}, ';'); end;
        reactomeID = getMetIDsInTargetLanguageFromInputLanguage(id,'MNX','reactome', metMNXxrefSourceID, metMNXxrefSourceMet, metMNXIDs);
        if ~isempty(reactomeID{1}) && ((~isempty(modelOut.metReactomeID{posMetInModel_i}) && overwrite) || isempty(modelOut.metReactomeID{posMetInModel_i})); modelOut.metReactomeID{posMetInModel_i} = strjoin(reactomeID{1}, ';'); end;
        rheaID = getMetIDsInTargetLanguageFromInputLanguage(id,'MNX','rheaG', metMNXxrefSourceID, metMNXxrefSourceMet, metMNXIDs);
        if ~isempty(rheaID{1}) && ((~isempty(modelOut.metRHEAID{posMetInModel_i}) && overwrite) || isempty(modelOut.metRHEAID{posMetInModel_i})); modelOut.metRHEAID{posMetInModel_i} = strjoin(rheaID{1}, ';'); end;
        sabiorkMID = getMetIDsInTargetLanguageFromInputLanguage(id,'MNX','sabiorkM', metMNXxrefSourceID, metMNXxrefSourceMet, metMNXIDs);
        if ~isempty(sabiorkMID{1}) && ((~isempty(modelOut.metSabiorkID{posMetInModel_i}) && overwrite) || isempty(modelOut.metSabiorkID{posMetInModel_i})); modelOut.metSabiorkID{posMetInModel_i} = strjoin(sabiorkMID{1}, ';'); end;
        slmID = getMetIDsInTargetLanguageFromInputLanguage(id,'MNX','slm', metMNXxrefSourceID, metMNXxrefSourceMet, metMNXIDs);
        if ~isempty(slmID{1}) && ((~isempty(modelOut.metSLMID{posMetInModel_i}) && overwrite) || isempty(modelOut.metSLMID{posMetInModel_i})); modelOut.metSLMID{posMetInModel_i} = strjoin(slmID{1}, ';'); end;
        
    end
end

for i = 1:length(posRxns)
    posRxnInModel_i = posRxns(i);
    if ~isempty(modelOut.rxnMetaNetXID{posRxnInModel_i})
        id = modelOut.rxnMetaNetXID(posRxnInModel_i);
        biggID = getRxnIDsInTargetLanguageFromInputLanguage_4_0(id,'MNX','biggR', rxnMNXxrefSourceID, rxnMNXxrefSourceRxn, rxnMNXxrefID);
        if ~isempty(biggID{1}) && ((~isempty(modelOut.rxnBiGGID{posRxnInModel_i}) && overwrite) || isempty(modelOut.rxnBiGGID{posRxnInModel_i})); modelOut.rxnBiGGID{posRxnInModel_i} = strjoin(biggID{1}, ';'); end;
        keggID =  getRxnIDsInTargetLanguageFromInputLanguage_4_0(id,'MNX','keggR', rxnMNXxrefSourceID, rxnMNXxrefSourceRxn, rxnMNXxrefID);
        if ~isempty(keggID{1}) && ((~isempty(modelOut.rxnKEGGID{posRxnInModel_i}) && overwrite) || isempty(modelOut.rxnKEGGID{posRxnInModel_i})); modelOut.rxnKEGGID{posRxnInModel_i} = strjoin(keggID{1}, ';'); end;
        seedID = getRxnIDsInTargetLanguageFromInputLanguage_4_0(id,'MNX','seedR', rxnMNXxrefSourceID, rxnMNXxrefSourceRxn, rxnMNXxrefID);
        if ~isempty(seedID{1}) && ((~isempty(modelOut.rxnSEEDID{posRxnInModel_i}) && overwrite) || isempty(modelOut.rxnSEEDID{posRxnInModel_i})); modelOut.rxnSEEDID{posRxnInModel_i} = strjoin(seedID{1}, ';'); end;
        metacycID = getRxnIDsInTargetLanguageFromInputLanguage_4_0(id,'MNX','metacycR', rxnMNXxrefSourceID, rxnMNXxrefSourceRxn, rxnMNXxrefID);
        if ~isempty(metacycID{1}) && ((~isempty(modelOut.rxnMetaCycID{posRxnInModel_i}) && overwrite) || isempty(modelOut.rxnMetaCycID{posRxnInModel_i})); modelOut.rxnMetaCycID{posRxnInModel_i} = strjoin(metacycID{1}, ';'); end;
        rheaID = getRxnIDsInTargetLanguageFromInputLanguage(id,'MNX','rheaR', rxnMNXxrefSourceID, rxnMNXxrefSourceRxn, rxnMNXIDs);
        if ~isempty(rheaID{1}) && ((~isempty(modelOut.rxnRHEAID{posRxnInModel_i}) && overwrite) || isempty(modelOut.rxnRHEAID{posRxnInModel_i})); modelOut.rxnRHEAID{posRxnInModel_i} = strjoin(rheaID{1}, ';'); end;
        sabiorkID = getRxnIDsInTargetLanguageFromInputLanguage(id,'MNX','sabiorkR', rxnMNXxrefSourceID, rxnMNXxrefSourceRxn, rxnMNXIDs);
        if ~isempty(sabiorkID{1}) && ((~isempty(modelOut.rxnSabiorkID{posRxnInModel_i}) && overwrite) || isempty(modelOut.rxnSabiorkID{posRxnInModel_i})); modelOut.rxnSabiorkID{posRxnInModel_i} = strjoin(sabiorkID{1}, ';'); end;
    end
end

end