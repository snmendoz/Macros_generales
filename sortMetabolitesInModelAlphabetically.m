function newModel = sortMetabolitesInModelAlphabetically(model)

newModel = model;

new_met_list = sort(model.mets);

i = zeros(10000000,1);
j = zeros(10000000,1);
s = zeros(10000000,1);

previous_count = 0;

for k = 1:length(new_met_list)
    disp(k)
    if k==582
        disp(k)
    end
    pos_k = find(strcmp(model.mets, new_met_list{k}));
    pos_rxns_k = find(model.S(pos_k,:));
    coefs = full(model.S(pos_k,pos_rxns_k));
    count = length(pos_rxns_k);
    i(previous_count+1:previous_count+count) = k*ones(count,1);
    j(previous_count+1:previous_count+count) = pos_rxns_k';
    s(previous_count+1:previous_count+count) = coefs';
    previous_count = previous_count+count;
end

i = i(1:previous_count);
j = j(1:previous_count);
s = s(1:previous_count);
new_S = sparse(i,j,s,length(model.mets),length(model.rxns));

newModel.S = new_S;
newModel.mets = new_met_list;
newModel.metNames = new_met_list;

end