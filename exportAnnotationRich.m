function exportAnnotationRich(model, modelFileName)

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


CBT_DB_FIELD_PROPS = [fields; additionalFields];
writeCbModel(model, 'fileName', modelFileName, 'format', 'sbml');
CBT_DB_FIELD_PROPS = fields;
end