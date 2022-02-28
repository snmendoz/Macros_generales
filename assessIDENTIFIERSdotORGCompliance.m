function [metAssessment, rxnAssessment, metFields, rxnFields] = assessIDENTIFIERSdotORGCompliance(model)

global CBT_DB_FIELD_PROPS

if isempty(CBT_DB_FIELD_PROPS)
    fileName = which('COBRA_structure_fields.csv');
    [raw] = descFileRead(fileName);
    
    fields = fieldnames(raw);
    for i = 1:numel(fields)
        %Convert everything to strings.
        raw.(fields{i}) = cellstr(raw.(fields{i}));
    end
    %Get the indices for database, qualifier and reference.
    relrows = cellfun(@(x) ischar(x) && ~isempty(x),raw.databaseid);
    relarray = [raw.databaseid(relrows),raw.qualifier(relrows),raw.Model_Field(relrows), raw.referenced_Field(relrows),raw.DBPatterns(relrows),raw.qualifierType(relrows)];
    dbInfo = cell(0,6);
    for i = 1:size(relarray)
        fieldRef = relarray{i,4}(1:end-1);
        dbs = strsplit(relarray{i,1},';');
        for db = 1:length(dbs)
            quals = strsplit(relarray{i,2},';');
            for qual = 1:length(quals)
                dbInfo(end+1,:) = {dbs{db},quals{qual},relarray{i,3},fieldRef,relarray{i,5},relarray{i,6}};
            end
        end
    end
    CBT_DB_FIELD_PROPS = dbInfo;
end
if size(CBT_DB_FIELD_PROPS,1)>18;
    CBT_DB_FIELD_PROPS = CBT_DB_FIELD_PROPS(1:18,:);
end
fields = CBT_DB_FIELD_PROPS;

% webread('https://registry.api.identifiers.org/restApi/namespaces/search/findByPrefix?prefix=inchikey')

additionalFields = {'metacyc.compound','is','metMetaCycID','met','.*','bioQualifier';...
    'metacyc.reaction','is','rxnMetaCycID','rxn','.*','bioQualifier';...
    'seed.compound','is','metSEEDID', 'met','^cpd\d+$','bioQualifier';...
    'seed.reaction','is','rxnSEEDID', 'rxn','^rxn\d+$','bioQualifier';...
    'kegg.glycan','is','metKEGGID', 'met','^G\d+$','bioQualifier';...
    'bigg.metabolite','is','metBiGGID','met','^[a-z_A-Z0-9]+$','bioQualifier';...
    'bigg.reaction','is','rxnBiGGID','rxn','^[a-z_A-Z0-9]+$','bioQualifier';...
    'rhea','is','rxnRHEAID','rxn','\d{5}$','bioQualifier';...
    'lipidmaps','is','metLipidMapsID','met','^LM(FA|GL|GP|SP|ST|PR|SL|PK)[0-9]{4}([0-9a-zA-Z]{4,6})?$','bioQualifier';...
    'inchikey','is','metInChIKey', 'met','^[A-Z]{14}\-[A-Z]{10}(\-[A-Z])?','bioQualifier';...
    'enviPath','is','metEnvipathID','met','^[\w^_]{8}-[\w^_]{4}-[\w^_]{4}-[\w^_]{4}-[\w^_]{12}\/[\w-]+\/[\w^_]{8}-[\w^_]{4}-[\w^_]{4}-[\w^_]{4}-[\w^_]{12}$','bioQualifier';...
    'sabiork.compound','is','metSabiorkID','met','^\d+$','bioQualifier';...
    'sabiork.reaction','is','rxnSabiorkID','rxn','^\d+$','bioQualifier';...
    'slm','is','metSLMID','met','^SLM:\d+$','bioQualifier'};

properties = [fields; additionalFields];
metFields = properties(startsWith(properties(:,3),'met'),3);
rxnFields = properties(startsWith(properties(:,3),'rxn'),3);
metExpressions = properties(startsWith(properties(:,3),'met'),5);
rxnExpressions = properties(startsWith(properties(:,3),'rxn'),5);

n_metFields = length(metFields);
n_rxnFields = length(rxnFields);
metAssessment = zeros(length(model.mets),n_metFields);
rxnAssessment = zeros(length(model.rxns),n_rxnFields);

for i = 1:n_metFields
    if isfield(model, metFields{i})
        infoField = eval(['model.' metFields{i}]);
        metAssessment(:,i) = cellfun(@(x) complyWithExpression(x, metExpressions{i}), infoField);
    end
end

for i = 1:n_rxnFields
    if isfield(model, rxnFields{i})
        infoField = eval(['model.' rxnFields{i}]);
        rxnAssessment(:,i) = cellfun(@(x) complyWithExpression(x, rxnExpressions{i}), infoField);
    end
end

end