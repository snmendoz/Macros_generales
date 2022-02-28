function [Presence_rxns,Reversibility, Presence_mets]=CreatePresenceAndReversibility(database)

load(database)

if ~exist('models.mat')
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
    end
    save(['modelNames_' num2str(length(modelNames)) '.mat'],'modelNames');
else
    load('models.mat')
end

Presence_mets=zeros(length(bigg.mets),length(models));
Presence_rxns=zeros(length(bigg.rxns),length(models));
Irreversible=NaN(length(bigg.rxns),length(models));
for i=1:length(models)
	disp(i)
    model_i=models{i};
    model_i_lb=model_i.lb;
    model_i_ub=model_i.ub;
    rxns_i=model_i.rxns;
    [~,ind_a,ind_b]=intersect(bigg.rxns,rxns_i);
    rev_pos=intersect(find(model_i_lb<0),find(model_i_ub>0));
    irev_izq_pos=intersect(find(model_i_lb<0),find(model_i_ub<=0));
    irev_der_pos=intersect(find(model_i_lb>=0),find(model_i_ub>0));
    
    [~,ind_1,~]=intersect(bigg.rxns,rxns_i(rev_pos));
    [~,ind_2,~]=intersect(bigg.rxns,rxns_i(irev_izq_pos));
    [~,ind_3,~]=intersect(bigg.rxns,rxns_i(irev_der_pos));
    
    if ~isempty(ind_2)
       disp('') 
    end
    
    Presence_rxns(ind_a,i)=ones(length(ind_a),1);
    Irreversible(ind_1,i)=0;
    Irreversible(ind_2,i)=-1;
    Irreversible(ind_3,i)=1;
    
    mets_i=model_i.mets;
    [~,ind_a,~]=intersect(bigg.mets,mets_i);
    Presence_mets(ind_a,i)=ones(length(ind_a),1);
end

save(['Presence_rxns_' database],'Presence_rxns');
save(['Irreversible_' database],'Irreversible');
save(['Presence_mets_' database],'Presence_mets');

end