function addSynonymsToBIGG(database)
cd('D:\Dropbox\Databases\BIGG')
load(database)
load('D:\Dropbox\Databases\MNX\rxnMNXIDs');
load('D:\Dropbox\Databases\MNX\rxnOtherIDs');

metHMD = cell(size(bigg.mets));
metMetaNetX = cell(size(bigg.mets));
metKEGG = cell(size(bigg.mets));
metSEED = cell(size(bigg.mets));
metBioCyc = cell(size(bigg.mets));
metReactome = cell(size(bigg.mets));
metLipidMaps = cell(size(bigg.mets));
metCHEBI = cell(size(bigg.mets));

rxnBioCyc = cell(size(bigg.rxns));
rxnKEGG = cell(size(bigg.rxns));
rxnSEED = cell(size(bigg.rxns));
rxnMetaNetX = cell(size(bigg.rxns));
rxnEC = cell(size(bigg.rxns));
rxnRHEA = cell(size(bigg.rxns));

posRevisar = zeros(1000,1);
cont_posRevisar = 0;
mets = regexprep(bigg.mets,'(.*)_(.*)$' ,'$1');
rxns = bigg.rxns;

for i = 1:length(bigg.mets)
    disp(i)
    [s, r] = system(['curl http://bigg.ucsd.edu/api/v2/universal/metabolites/' mets{i}]);
    if s==0
        posF = strfind(r, '/kegg.compound');
        kegg = {};
        if ~isempty(posF)
            for j = 1:length(posF)
                posSlash = strfind(r, '/');
                posComma = strfind(r, '"');
                posComma1 = posComma(posComma>posF(j));
                posSlash1 = posSlash(posSlash>posF(j));
                r1 = r(posSlash1(1)+1:posComma1(1)-1);
                kegg = [kegg, r1];
            end
        end
        metKEGG{i} = strjoin(unique(kegg),',');
        
        posF = strfind(r, '/seed.compound');
        seed = {};
        if ~isempty(posF)
            for j = 1:length(posF)
                posSlash = strfind(r, '/');
                posComma = strfind(r, '"');
                posComma1 = posComma(posComma>posF(j));
                posSlash1 = posSlash(posSlash>posF(j));
                r1 = r(posSlash1(1)+1:posComma1(1)-1);
                seed = [seed, r1];
            end
        end
        metSEED{i} = strjoin(unique(seed),',');
        
        posF = strfind(r, '/metanetx.chemical');
        mnx = {};
        if ~isempty(posF)
            for j = 1:length(posF)
                posSlash = strfind(r, '/');
                posComma = strfind(r, '"');
                posComma1 = posComma(posComma>posF(j));
                posSlash1 = posSlash(posSlash>posF(j));
                r1 = r(posSlash1(1)+1:posComma1(1)-1);
                mnx = [mnx, r1];
            end
        end
        metMetaNetX{i} = strjoin(unique(mnx),',');
        
        posF = strfind(r, '/biocyc');
        biocyc = {};
        if ~isempty(posF)
            for j = 1:length(posF)
                posSlash = strfind(r, '/');
                posComma = strfind(r, '"');
                posComma1 = posComma(posComma>posF(j));
                posSlash1 = posSlash(posSlash>posF(j));
                r1 = r(posSlash1(1)+1:posComma1(1)-1);
                biocyc = [biocyc, regexprep(r1,'META:','')];
            end
        end
        metBioCyc{i} = strjoin(unique(biocyc),',');
        
        posF = strfind(r, '/hmdb');
        hmdb = {};
        if ~isempty(posF)
            for j = 1:length(posF)
                posSlash = strfind(r, '/');
                posComma = strfind(r, '"');
                posComma1 = posComma(posComma>posF(j));
                posSlash1 = posSlash(posSlash>posF(j));
                r1 = r(posSlash1(1)+1:posComma1(1)-1);
                hmdb = [hmdb, r1];
            end
        end
        metHMD{i} = strjoin(unique(hmdb),',');
        
        posF = strfind(r, '/lipidmaps');
        lipidmaps = {};
        if ~isempty(posF)
            for j = 1:length(posF)
                posSlash = strfind(r, '/');
                posComma = strfind(r, '"');
                posComma1 = posComma(posComma>posF(j));
                posSlash1 = posSlash(posSlash>posF(j));
                r1 = r(posSlash1(1)+1:posComma1(1)-1);
                lipidmaps = [lipidmaps, r1];
            end
        end
        metLipidMaps{i} = strjoin(unique(lipidmaps),',');
        
        posF = strfind(r, '/chebi');
        chebi = {};
        if ~isempty(posF)
            for j = 1:length(posF)
                posSlash = strfind(r, '/');
                posComma = strfind(r, '"');
                posComma1 = posComma(posComma>posF(j));
                posSlash1 = posSlash(posSlash>posF(j));
                r1 = r(posSlash1(1)+1:posComma1(1)-1);
                chebi = [chebi, regexprep(r1,'CHEBI:','')];
            end
        end
        metCHEBI{i} = strjoin(unique(chebi),',');
        
        posF = strfind(r, '/reactome');
        reactome = {};
        if ~isempty(posF)
            for j = 1:length(posF)
                posSlash = strfind(r, '/');
                posComma = strfind(r, '"');
                posComma1 = posComma(posComma>posF(j));
                posSlash1 = posSlash(posSlash>posF(j));
                r1 = r(posSlash1(1)+1:posComma1(1)-1);
                reactome = [reactome, r1];
            end
        end
        metReactome{i} = strjoin(unique(reactome),',');
        
    else
        disp('')
        cont_posRevisar = cont_posRevisar + 1;
        posRevisar(cont_posRevisar) = i;
    end
end

for i = 1:length(bigg.rxns)
    [s, r] = system(['curl http://bigg.ucsd.edu/api/v2/universal/reactions/' rxns{i}]);
    disp(i)
    posF = strfind(r, '/biocyc');
    biocyc = {};
    if ~isempty(posF)
        for j = 1:length(posF)
            posSlash = strfind(r, '/');
            posComma = strfind(r, '"');
            posComma1 = posComma(posComma>posF(j));
            posSlash1 = posSlash(posSlash>posF(j));
            r1 = r(posSlash1(1)+1:posComma1(1)-1);
            biocyc = [biocyc, regexprep(r1,'META:','')];
        end
    end
    rxnBioCyc{i} = strjoin(unique(biocyc),',');
    
    posF = strfind(r, '/kegg.reaction');
    kegg = {};
    if ~isempty(posF)
        for j = 1:length(posF)
            posSlash = strfind(r, '/');
            posComma = strfind(r, '"');
            posComma1 = posComma(posComma>posF(j));
            posSlash1 = posSlash(posSlash>posF(j));
            r1 = r(posSlash1(1)+1:posComma1(1)-1);
            kegg = [kegg, r1];
        end
    end
    rxnKEGG{i} = strjoin(unique(kegg),',');
    
    %         posF = strfind(r, '/seed.reaction');
    %         seed = {};
    %         if ~isempty(posF)
    %             for j = 1:length(posF)
    %                 posSlash = strfind(r, '/');
    %                 posComma = strfind(r, '"');
    %                 posComma1 = posComma(posComma>posF(j));
    %                 posSlash1 = posSlash(posSlash>posF(j));
    %                 r1 = r(posSlash1(1)+1:posComma1(1)-1);
    %                 seed = [seed, r1];
    %             end
    %         end
    %         rxnSEED{i} = strjoin(unique(seed),',');
    
    posF = strfind(r, '/metanetx.reaction');
    mnx = {};
    if ~isempty(posF)
        for j = 1:length(posF)
            posSlash = strfind(r, '/');
            posComma = strfind(r, '"');
            posComma1 = posComma(posComma>posF(j));
            posSlash1 = posSlash(posSlash>posF(j));
            r1 = r(posSlash1(1)+1:posComma1(1)-1);
            mnx = [mnx, r1];
        end
    end
    rxnMetaNetX{i} = strjoin(unique(mnx),',');
    
    if ~isempty(rxnMetaNetX{i})
        id = getRxnIDsInTargetLanguageFromInputLanguage(rxnMetaNetX(i), {'MNX'},{'seed'},rxnOtherIDs, rxnMNXIDs);
        id = id{1};
        if ~isempty(id)
            rxnSEED{i} = id{1};
        end
    end
    
    posF = strfind(r, '/ec');
    ec = {};
    if ~isempty(posF)
        for j = 1:length(posF)
            posSlash = strfind(r, '/');
            posComma = strfind(r, '"');
            posComma1 = posComma(posComma>posF(j));
            posSlash1 = posSlash(posSlash>posF(j));
            r1 = r(posSlash1(1)+1:posComma1(1)-1);
            ec = [ec, r1];
        end
    end
    rxnEC{i} = strjoin(unique(ec),',');
    
    posF = strfind(r, '/rhea');
    rhea = {};
    if ~isempty(posF)
        for j = 1:length(posF)
            posSlash = strfind(r, '/');
            posComma = strfind(r, '"');
            posComma1 = posComma(posComma>posF(j));
            posSlash1 = posSlash(posSlash>posF(j));
            r1 = r(posSlash1(1)+1:posComma1(1)-1);
            rhea = [rhea, r1];
        end
    end
    rxnRHEA{i} = strjoin(unique(rhea),',');
    
end

posRevisar = posRevisar(1:cont_posRevisar);
save('posRevisar', 'posRevisar')

bigg.metKEGGID = metKEGG;
bigg.metSEEDID = metSEED;
bigg.metMetaNetXID = metMetaNetX;
bigg.metMetaCycID = metBioCyc;
bigg.metReactomeID = metReactome;
bigg.metLipidMapsID = metLipidMaps;
bigg.metHMDID = metHMD;
bigg.metChEBIID = metCHEBI;
bigg.metBiGGID = removeCompartmentFromMets(bigg.mets);

bigg.rxnBiGGID = bigg.rxns;
bigg.rxnMetaCycID = rxnBioCyc;
bigg.rxnKEGGID = rxnKEGG;
bigg.rxnSEEDID = rxnSEED;
bigg.rxnMetaNetXID = rxnMetaNetX;
bigg.rxnECNumbers = rxnEC;
bigg.rxnRHEAID = rxnRHEA;


for i = 1:length(bigg.metChEBIID); if isempty(bigg.metChEBIID{i}); bigg.metChEBIID{i} =''; end; end;
for i = 1:length(bigg.metKEGGID); if isempty(bigg.metKEGGID{i}); bigg.metKEGGID{i} =''; end; end;
for i = 1:length(bigg.metSEEDID); if isempty(bigg.metSEEDID{i}); bigg.metSEEDID{i} =''; end; end;
for i = 1:length(bigg.metMetaNetXID); if isempty(bigg.metMetaNetXID{i}); bigg.metMetaNetXID{i} =''; end; end;
for i = 1:length(bigg.metMetaCycID); if isempty(bigg.metBioCycID{i}); bigg.metBioCycID{i} =''; end; end;
for i = 1:length(bigg.metHMDID); if isempty(bigg.metHMDID{i}); bigg.metHMDID{i} =''; end; end;
for i = 1:length(bigg.metReactomeID); if isempty(bigg.metReactomeID{i}); bigg.metReactomeID{i} =''; end; end;
for i = 1:length(bigg.metLipidMapsID); if isempty(bigg.metLipidMapsID{i}); bigg.metLipidMapsID{i} =''; end; end;

for i = 1:length(bigg.rxnMetaCycID); if isempty(bigg.rxnBioCycID{i}); bigg.rxnBioCycID{i} =''; end; end;
for i = 1:length(bigg.rxnKEGGID); if isempty(bigg.rxnKEGGID{i}); bigg.rxnKEGGID{i} =''; end; end;
for i = 1:length(bigg.rxnSEEDID); if isempty(bigg.rxnSEEDID{i}); bigg.rxnSEEDID{i} =''; end; end;
for i = 1:length(bigg.rxnMetaNetXID); if isempty(bigg.rxnMetaNetXID{i}); bigg.rxnMetaNetXID{i} =''; end; end;
for i = 1:length(bigg.rxnECNumbers); if isempty(bigg.rxnECID{i}); bigg.rxnECID{i} =''; end; end;
for i = 1:length(bigg.rxnRHEAID); if isempty(bigg.rxnRHEAID{i}); bigg.rxnRHEAID{i} =''; end; end;

for i = 1:length(bigg.metChEBIID)
    if ~isempty(bigg.metChEBIID{i})
        bigg.metChEBIID{i} = ['CHEBI:' regexprep(bigg.metChEBIID{i},',',';CHEBI:')];
    end
end
save(database, 'bigg');