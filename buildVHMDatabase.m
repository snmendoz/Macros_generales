function buildVHMDatabase

load('modelsID')
models = cell(773, 1);
cd('C:\Users\notebook\Desktop\Trabajo\Gut\AllModels\MAT')
if exist('C:\Users\notebook\Desktop\Trabajo\Gut\AllModels\MAT\models.mat', 'file')~=2       
    for i=1:length(modelsID)
        try
            load(modelsID{i});
        catch
            disp(['https://webdav-r3lab.uni.lu/public/msp/AGORA/mat/' modelsID{i} '.mat'])
            load(modelsID{i});
        end
        models{i}=model;
        clear model; 
    end
    save('models', 'models','-v7.3')
else
    load('models2.mat')
end

for i=1:length(models)
    model_i=models{i};
    pos_biomass1=find(cellfun(@isempty,strfind(model_i.rxns,'biomass'))==0);
    pos_biomass2=find(cellfun(@isempty,strfind(model_i.rxns,'Biomass'))==0);
    pos_biomass3=find(cellfun(@isempty,strfind(model_i.rxns,'BIOMASS'))==0);
    pos_biomass=[pos_biomass1;pos_biomass2;pos_biomass3];
    if isempty(pos_biomass)
        disp('')
    else
        model_i=removeRxns(model_i,model_i.rxns(pos_biomass));
    end    
    
    if i==1;
        vhm=model_i;
        vhm.description='vhm';
    else
        rxns_i=model_i.rxns;
        [Faltantes,ind]=setdiff(rxns_i,vhm.rxns);
        ec=EscribirRxn_cobraFormat(model_i,ind);
        [pos_mets,~]=find(model_i.S(:,ind));
        pos_mets=unique(pos_mets);
        MetFaltantes=setdiff(model_i.mets(pos_mets),vhm.mets);
        [~,~,ind_m]=intersect(MetFaltantes,model_i.mets);
        for j=1:length(MetFaltantes)
            vhm=addMetabolite(vhm,MetFaltantes(j),'metName', model_i.metNames(ind_m(j)));
        end
        for j=1:length(ec)
%             fprintf('i: %.0f, j: %.0f\n',i,j)
                if ReactionInModel(vhm, ec{j})
                    [~, pos] = ReactionInModel(vhm, ec{j});
                    disp('repeatedIDs')
                else
                    vhm=addReaction(vhm,Faltantes{j},ec{j});
                end
            
        end
    end
end
save('vhm','vhm');

end