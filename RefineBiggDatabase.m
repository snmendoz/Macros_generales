function names = RefineBiggDatabase(database)
load(database)
comparments=cell(length(bigg.mets),1);
metsToRemove={};
for i=1:length(bigg.mets)
    if ~isempty(strfind(bigg.mets{i},'['))
        comparments{i}=bigg.mets{i}(end-1);
    else
        if strcmp(bigg.mets{i}(end-1),'_')
            comparments{i}=bigg.mets{i}(end);
            aux=[bigg.mets{i}(1:end-2) '[' bigg.mets{i}(end) ']'];
            if ~isempty(strmatch(aux,bigg.mets,'exact'))
                %
                pos_rxns=find(bigg.S(i,:));
                posMet=strmatch(aux,bigg.mets,'exact');
                bigg.S(posMet,pos_rxns)=bigg.S(i,pos_rxns);
                bigg.S(i,pos_rxns)=0;
                metsToRemove=union(metsToRemove,bigg.mets(i));
            else
                bigg.mets{i}=aux;
            end
        else 
            comparments{i}=bigg.mets{i}(end-1:end);
            bigg.mets{i}=[bigg.mets{i}(1:end-3) '[' bigg.mets{i}(end-1:end) ']'];
        end
    end
end
bigg=removeMetabolites(bigg,metsToRemove);
% listComparments=unique(comparments);
% for i=1:length(listComparments)
%     bigg.mets=regexprep(bigg.mets,['_' listComparments{i}],['[' listComparments{i} ']']);
% end

%fix reversibility of reactions
load(['Presence_rxns_' database]);
load(['Irreversible_' database]);

for i=1:length(bigg.rxns)
    pos_rxns=find(Presence_rxns(i,:));
    rate_ir_izq=100*(length(find(Irreversible(i,pos_rxns)==-1))/length(pos_rxns));
    rate_ir_der=100*(length(find(Irreversible(i,pos_rxns)==1))/length(pos_rxns));
    if rate_ir_izq==100
        bigg.lb(i)=-1000;
        bigg.ub(i)=0;
    elseif rate_ir_der==100
        bigg.lb(i)=0;
        bigg.ub(i)=1000;
    else
        bigg.lb(i)=-1000;
        bigg.ub(i)=1000;
    end
    
end

bigg=rmfield(bigg,'genes');
bigg=rmfield(bigg,'grRules');
bigg=rmfield(bigg,'subSystems');
bigg=rmfield(bigg,'rxnGeneMat');

save([database '_fixed'],'bigg');
posRxnsWithOneCompound=find(sum(bigg.S~=0,1)==1);
bigg=removeRxns(bigg,bigg.rxns(posRxnsWithOneCompound));
save([database '_fixed_forMeneco'],'bigg');

names{1} = [database '_fixed'];
names{2} = [database '_fixed_forMeneco'];

% current=pwd;
% cd('C:\Users\notebook\Desktop\Trabajo\cobra_antiguo\io')
% writeCbModel2(bigg,'sbml',[database '_fixed'])
% writeCbModel2(biggWithoutExDmSkRxns,'sbml',[database '_fixed_Meneco'])
% cd(current)
end