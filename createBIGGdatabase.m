function names = createBIGGdatabase
cd('D:\Dropbox\Databases\BIGG')

if exist('models.mat','file')~=2 || exist('models_refined.mat','file')~=2
    
    fi=fopen('modelNames.txt');
    tline = fgetl(fi);
    modelNames=[];
    while ~isempty(tline)
        modelNames=[modelNames;{tline}];
        tline = fgetl(fi);
        if tline==-1
            break;
        end
    end
    fclose(fi);
    
    for i=1:length(modelNames)
        load(modelNames{i});
        models{i}=eval(modelNames{i});
        models_refined{i} = refineModel(models{i});

    end
    save('models','models','-v7.3')
    save('models_refined','models_refined','-v7.3')
else
    load('models.mat')
    load('models_refined.mat')
    
    fi=fopen('modelNames.txt');
    tline = fgetl(fi);
    modelNames=[];
    while ~isempty(tline)
        modelNames=[modelNames;{tline}];
        tline = fgetl(fi);
        if tline==-1
            break;
        end
    end
    fclose(fi);
    n_models = length(modelNames);
    if length(models) ~= n_models
        for i=1:length(modelNames)
            load(modelNames{i});
            models{i}=eval(modelNames{i});
            models_refined{i} = refineModel(models{i});
            
        end
        save('models','models','-v7.3')
        save('models_refined','models_refined','-v7.3')
    end
end

if exist(['bigg_' num2str(length(models)) '.mat'],'file')==2
    load(['bigg_' num2str(length(models))])
else
    
    for i=1:length(models)
        model_i=models{1,i};
        pos_biomass1=find(cellfun(@isempty,strfind(model_i.rxns,'biomass'))==0);
        pos_biomass2=find(cellfun(@isempty,strfind(model_i.rxns,'Biomass'))==0);
        pos_biomass3=find(cellfun(@isempty,strfind(model_i.rxns,'BIOMASS'))==0);
        pos_biomass=[pos_biomass1;pos_biomass2;pos_biomass3];
        if isempty(pos_biomass)
            disp('')
        else
            model_i=removeRxns(model_i,model_i.rxns(pos_biomass));
        end
        model_i.mets = regexprep(model_i.mets,{'\[intermembrane space of mitochondria\]','\[thylakoid membrane\]','\[cytosolic membrane\]','\[carboxyzome\]','\[mitochondrial intermembrane\]'},{'','','','',''});
        
        if i==1;
            bigg=model_i;
            bigg.description='bigg';
        else
            rxns_i=model_i.rxns;
            [Faltantes,ind]=setdiff(rxns_i,bigg.rxns);
            if ~isempty(Faltantes)
                ec=getRxn_cobraFormat(model_i,ind);
                [pos_mets,~]=find(model_i.S(:,ind));
                pos_mets=unique(pos_mets);
                MetFaltantes=setdiff(model_i.mets(pos_mets),bigg.mets);
                [~,~,ind_m]=intersect(MetFaltantes,model_i.mets);
                for j=1:length(MetFaltantes)
                    bigg=addMetabolite(bigg,MetFaltantes{j},model_i.metNames{ind_m(j)});
                end
                for j=1:length(ec)
                    if i ==18 && j==1
                        disp('')
                    end
                    fprintf('i: %.0f, j: %.0f\n',i,j)
                    bigg=addReaction(bigg,Faltantes{j},ec{j});
                end
            end
        end
    end
    
    save(['bigg_' num2str(length(models))],'bigg');
end

if exist(['bigg_' num2str(length(models)) '_Lp_Af.mat'],'file')==2

else
    
    %Agregar L.plantarum
    load('D:\Dropbox\Databases\SystemBioinformaticsModels\Lactobacillus plantarum\Lactobacillus_plantarum_WCFS1_19_Jul_2018_19_22_52.mat');
    model_i=model;
    
    pos_biomass1=find(cellfun(@isempty,strfind(model_i.rxns,'biomass'))==0);
    pos_biomass2=find(cellfun(@isempty,strfind(model_i.rxns,'Biomass'))==0);
    pos_biomass3=find(cellfun(@isempty,strfind(model_i.rxns,'BIOMASS'))==0);
    pos_biomass=[pos_biomass1;pos_biomass2;pos_biomass3];
    if isempty(pos_biomass)
        disp('')
    else
        model_i=removeRxns(model_i,model_i.rxns(pos_biomass));
    end
    
    rxns_i=model_i.rxns;
    [Faltantes,ind]=setdiff(rxns_i,bigg.rxns);
    ec=getRxn_cobraFormat(model_i,ind);
    [pos_mets,~]=find(model_i.S(:,ind));
    pos_mets=unique(pos_mets);
    MetFaltantes=setdiff(model_i.mets(pos_mets),bigg.mets);
    [~,~,ind_m]=intersect(MetFaltantes,model_i.mets);
    for j=1:length(MetFaltantes)
        bigg=addMetabolite(bigg,MetFaltantes{j},model_i.metNames{ind_m(j)});
    end
    for j=1:length(ec)
        bigg=addReaction(bigg,Faltantes{j},ec{j});
    end
    
    %Agregar A.ferroxidans
    load('mmc3_bigg_compliant');
    model_i=mmc3_bigg_compliant;
    
    pos_biomass1=find(cellfun(@isempty,strfind(model_i.rxns,'biomass'))==0);
    pos_biomass2=find(cellfun(@isempty,strfind(model_i.rxns,'Biomass'))==0);
    pos_biomass3=find(cellfun(@isempty,strfind(model_i.rxns,'BIOMASS'))==0);
    pos_biomass=[pos_biomass1;pos_biomass2;pos_biomass3];
    if isempty(pos_biomass)
        disp('')
    else
        model_i=removeRxns(model_i,model_i.rxns(pos_biomass));
    end
    
    rxns_i=model_i.rxns;
    [Faltantes,ind]=setdiff(rxns_i,bigg.rxns);
    ec=getRxn_cobraFormat(model_i,ind);
    [pos_mets,~]=find(model_i.S(:,ind));
    pos_mets=unique(pos_mets);
    MetFaltantes=setdiff(model_i.mets(pos_mets),bigg.mets);
    [~,~,ind_m]=intersect(MetFaltantes,model_i.mets);
    for j=1:length(MetFaltantes)
        bigg=addMetabolite(bigg,MetFaltantes{j},model_i.metNames{ind_m(j)});
    end
    for j=1:length(ec)
        bigg=addReaction(bigg,Faltantes{j},ec{j});
    end
    save(['bigg_' num2str(length(models)) '_Lp_Af'],'bigg');
end

names{1} = ['bigg_' num2str(length(models))];
names{2} = ['bigg_' num2str(length(models)) '_Lp_Af'];

end

function model = refineModel(model)
if length(find(cellfun(@isempty,strfind(model.mets,'['))==0)) ~= length(model.mets)
    mets = model.mets;
    for i = 1:length(mets)
        if strcmp(mets{i}(end-1),'_')
            mets{i} = [mets{i}(1:end-2) '[' mets{i}(end) ']'];
        elseif ~isempty(strfind(mets{i},'_im[intermembrane space of mitochondria]')) || ~isempty(strfind(mets{i},'_cx[carboxyzome')) || ~isempty(strfind(mets{i},'_um[thylakoid membrane]')) || ~isempty(strfind(mets{i},'_cm[cytosolic membrane]')) || ~isempty(strfind(mets{i},'_mm[mitochondrial intermembrane]'))
            mets{i} = regexprep(mets{i}, '_(.*)\[.*\]$','[$1]');
%             mets{i} = regexprep(mets{i}, {'_im[intermembrane space of mitochondria]','_cx[carboxyzome','_um[thylakoid membrane]','_cm[cytosolic membrane]'}, {'[im]','[cx]','[um]','[cm]'});
        else
            disp('caso raro')
        end
    end
    model.mets = mets;
end
end