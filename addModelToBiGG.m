function addModelToBiGG(databaseName, newDatabaseName, modelPath)

load(['D:\Dropbox\Databases\BIGG\' databaseName '.mat'])
load([modelPath '.mat'])
bigg = bigg;
bigg.mets = regexprep(bigg.mets,'_cx\[carboxyzome\]','_cx');

pos_biomass1=find(cellfun(@isempty,strfind(model.rxns,'biomass'))==0);
pos_biomass2=find(cellfun(@isempty,strfind(model.rxns,'Biomass'))==0);
pos_biomass3=find(cellfun(@isempty,strfind(model.rxns,'BIOMASS'))==0);
pos_biomass=[pos_biomass1;pos_biomass2;pos_biomass3];
if isempty(pos_biomass)
    disp('')
else
    model=removeRxns(model,model.rxns(pos_biomass));
end

rxns_i=model.rxns;
[diff,ind]=setdiff(rxns_i,bigg.rxns);
ec=getRxn_cobraFormat(model,ind);
[pos_mets,~]=find(model.S(:,ind));
pos_mets=unique(pos_mets);
Metdiff=setdiff(model.mets(pos_mets),bigg.mets);
[~,~,ind_m]=intersect(Metdiff,model.mets);
for j=1:length(Metdiff)
    bigg=addMetabolite(bigg,Metdiff{j},model.metNames{ind_m(j)});
end
for j=1:length(ec)
    bigg=addReaction(bigg,diff{j},ec{j});
end
cd('D:\Dropbox\Databases\BIGG')
save(['D:\Dropbox\Databases\BIGG\' newDatabaseName '.mat'], 'bigg');

posRxnsWithOneCompound=find(sum(bigg.S~=0,1)==1);
bigg=removeRxns(bigg,bigg.rxns(posRxnsWithOneCompound));
save([newDatabaseName '_fixed_forMeneco'],'bigg');
exportDataBaseToMeneco('bigg_85_more_fixediNF517_fixed_forMeneco')


end