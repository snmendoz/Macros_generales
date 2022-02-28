function essentials = predictEssentialNutrients(model, threshold)


model.lb(find(cellfun(@isempty, strfind(model.rxns, 'EX_'))==0)) = -1*ones(size(find(cellfun(@isempty, strfind(model.rxns, 'EX_'))==0)));
model.ub(find(cellfun(@isempty, strfind(model.rxns, 'EX_'))==0)) = 1*ones(size(find(cellfun(@isempty, strfind(model.rxns, 'EX_'))==0)));


fba = optimizeCbModel(model);
if nargin < 2;
   threshold = fba.f*0.1; 
end
posEx = find(cellfun(@isempty, strfind(model.rxns, 'EX_'))==0);
n_essentials = 0;
essentials = {};

for i = 1:length(posEx)
   disp(i)
   if i ==5
      disp('') 
   end
   model2 = changeRxnBounds(model,model.rxns(posEx(i)),0,'l') ;
   fba2 = optimizeCbModel(model2);
   if fba2.f < threshold
       n_essentials = n_essentials + 1;
       essentials = union(essentials, model.rxns(posEx(i)));
   end
end

end