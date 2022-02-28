function summary

ef =  readCbModel('efaecalis_final_checked.sbml');
model = ef;
%determine % of reactions with gene associations
pga = 100*(length(find(cellfun(@isempty,model.rules)==0))/length(model.rules));

%detemine if there is an artificial generation of ATP
% model2 = model; 
% model2 = changeRxnBounds(model2, model.rxns(strmatch('Exchange', model.rxns)), 0, 'b');

%determine number of bloqued reactions
blocked = zeros(size(model.rxns));
for i = 1:length(model.rxns)
   model2 =  changeObjective(model, model.rxns(i));
   fmin = optimizeCbModel(model2);
   fmax = optimizeCbModel(model2, 'min');
   if fmin.f<10^-9 && fmax.f<10^-9
       blocked(i) = 1;
   end
end
n_blocked = length(find(blocked));
pnb = 100*n_blocked/length(model.rxns);

%determine list of essential rxns and genes
essentials = zeros(size(model.rxns));
for i = 1:length(model.rxns)
   model2 =  changeRxnBounds(model, model.rxns(i), 0, 'b');
   fmax = optimizeCbModel(model2);
   if fmax.f<10^-9
       essentials(i) = 1;
   end
end
rxns_ess = model.rxns(find(essentials));


end