function [interaction, alone1, alone2, best1, best2] = DetermineInteractionUsingDual(runID, model1, model2, obj, medium, species)

%define medium for independent growth
model1 = changeRxnBounds(model1, model1.rxns(cellfun(@isempty, strfind(model1.rxns, 'EX_'))==0), 0, 'l');
model2 = changeRxnBounds(model2, model2.rxns(cellfun(@isempty, strfind(model2.rxns, 'EX_'))==0), 0, 'l');

model1 = changeRxnBounds(model1, medium.Rxns, -medium.Rates, 'l');
model2 = changeRxnBounds(model2, medium.Rxns, -medium.Rates, 'l');

%define independet growth rates
fba1 = optimizeCbModel(model1);
fba2 = optimizeCbModel(model2);

if fba1.f == 0 || fba2.f == 0
    warning('growth rate is zero')
    disp('')
end

%determine common compounds and independently consumed compounds
shared = findSharedExtracellularMetabolites(model1, model2);
shared = regexprep(shared,'\[e\]','');

onlyModel1 = setdiff(model1.mets(cellfun(@isempty,strfind(model1.mets, '[e]'))==0), model2.mets(cellfun(@isempty,strfind(model2.mets, '[e]'))==0));
onlyModel2 = setdiff(model2.mets(cellfun(@isempty,strfind(model2.mets, '[e]'))==0), model1.mets(cellfun(@isempty,strfind(model1.mets, '[e]'))==0));
indMets = {onlyModel1, onlyModel2};

lb_pool1 = zeros(size(onlyModel1));
for j = 1:length(onlyModel1)
    pos_met = find(strcmp(model1.mets, onlyModel1(j)));
    pos_ex = intersect(find(sum(model1.S~=0,1)==1),find(model1.S(pos_met,:)));
    if ismember(model1.rxns(pos_ex), medium.Rxns)
        [bol, pos] = ismember(model1.rxns(pos_ex), medium.Rxns);
        lb_pool1(j) = -medium.Rates(pos);
    end 
end

lb_pool2 = zeros(size(onlyModel2));
for j = 1:length(onlyModel2)
    pos_met = find(strcmp(model2.mets, onlyModel2(j)));
    pos_ex = intersect(find(sum(model2.S~=0,1)==1),find(model2.S(pos_met,:)));
    if ismember(model2.rxns(pos_ex), medium.Rxns)
        [bol, pos] = ismember(model2.rxns(pos_ex), medium.Rxns);
        lb_pool2(j) = -medium.Rates(pos);
    end 
end
 
lb = zeros(size(shared));
for j = 1:length(shared)
    met_i = [shared{j} '[e]'];
    pos_met1 = find(strcmp(model1.mets, met_i));
    pos_ex1 = intersect(find(sum(model1.S~=0,1)==1),find(model1.S(pos_met1,:)));
    
    pos_met2 = find(strcmp(model2.mets, met_i));
    pos_ex2 = intersect(find(sum(model2.S~=0,1)==1),find(model2.S(pos_met2,:)));
    
    if strcmp(model1.rxns(pos_ex1), model2.rxns(pos_ex2))
        if ismember(model1.rxns(pos_ex1), medium.Rxns)
            [bol, pos] = ismember(model1.rxns(pos_ex1), medium.Rxns);
            lb(j) = -medium.Rates(pos);
        end
    else
        if ismember(model1.rxns(pos_ex1), medium.Rxns)
            [bol, pos] = ismember(model1.rxns(pos_ex1), medium.Rxns);
            lb(j) = -medium.Rates(pos);
        end
        if ismember(model2.rxns(pos_ex2), medium.Rxns)
            [bol, pos] = ismember(model2.rxns(pos_ex2), medium.Rxns);
            lb(j) = -medium.Rates(pos);
        end
        disp('')
    end    
end

%merge Models
models = {model1, model2};
lb_pool = {lb_pool1, lb_pool2};
% lb = -1*ones(length(shared),1);
[mergedModel, exchangeRxns] = mergeModels(models, shared, 'lbPool', lb, 'indLb', lb_pool, 'indMets', indMets, 'dual', 1);

%define medium constraints
for i = 1:length(exchangeRxns)
    if ismember(regexprep(exchangeRxns{i}, {'model1_', 'model2_'}, {'',''}), medium.Rxns)
        [bol, pos] = ismember(regexprep(exchangeRxns{i}, {'model1_', 'model2_'}, {'',''}), medium.Rxns);
        mergedModel = changeRxnBounds(mergedModel, exchangeRxns{i}, -medium.Rates(pos), 'l');
    else
        mergedModel = changeRxnBounds(mergedModel, exchangeRxns{i}, 0, 'l');
    end
end

%verify biomass optimization

posObj = [strmatch(['model1_' obj{1}], mergedModel.rxns); strmatch(['model2_' obj{2}], mergedModel.rxns)];

mergedModel1 = changeObjective(mergedModel, ['model1_' obj{1}]);
fba1b = optimizeCbModel(mergedModel1);
alone1 = fba1b.f;
mergedModel2 = changeObjective(mergedModel, ['model2_' obj{2}]);
fba2b = optimizeCbModel(mergedModel2);
alone2 = fba2b.f;

if fba1b.f == 0 || fba2b.f == 0
    warning('growth rate is zero')
    disp('')
end

if fba1b.f < 0.99*fba1.f || fba1b.f > 1.01*fba1.f || fba2b.f < 0.99*fba2.f || fba2b.f > 1.01*fba2.f
    warning('different growth rates')
    disp('')
end

%%
% xs_2 = [];
% ys_2 = [];
% 
% mergedModelParetoManual.c(posObj(1)) = 1;
% mergedModelParetoManual.c(posObj(2)) = 0;
% fbamin = optimizeCbModel(mergedModelParetoManual, 'min');
% fbamax = optimizeCbModel(mergedModelParetoManual, 'max');
% 
% span = linspace(fbamin.f, fbamax.f, 100);
% mergedModelParetoManual.c(posObj(1)) = 0;
% mergedModelParetoManual.c(posObj(2)) = 1;
% 
% for i = span
%     mergedModelParetoManual2 = changeRxnBounds(mergedModelParetoManual, mergedModelParetoManual.rxns(posObj(1)), i, 'b');
%     fba_aux = optimizeCbModel(mergedModelParetoManual2);
%     xs_2 = [xs_2; fba_aux.f];
% end
% figure;
% plot(span,xs_2);
% 
% 
% mergedModelParetoManual.c(posObj(1)) = 0;
% mergedModelParetoManual.c(posObj(2)) = 1;
% fbamin = optimizeCbModel(mergedModelParetoManual, 'min');
% fbamax = optimizeCbModel(mergedModelParetoManual, 'max');
% 
% span = linspace(fbamin.f, fbamax.f, 100);
% mergedModelParetoManual.c(posObj(1)) = 1;
% mergedModelParetoManual.c(posObj(2)) = 0;
% 
% for i = span
%     mergedModelParetoManual2 = changeRxnBounds(mergedModelParetoManual, mergedModelParetoManual.rxns(posObj(2)), i, 'b');
%     fba_aux = optimizeCbModel(mergedModelParetoManual2);
%     ys_2 = [ys_2; fba_aux.f];
% end
% figure;
% plot(span,ys_2);

%%

%determine bestpoint
mergedModelTogether = mergedModel;
mergedModelTogether.c(strcmp(mergedModelTogether.rxns, ['model1_' obj{1}])) = 1;
mergedModelTogether.c(strcmp(mergedModelTogether.rxns, ['model2_' obj{2}])) = 1;

together = optimizeCbModel(mergedModelTogether);
best1 = together.x(find(strcmp(mergedModelTogether.rxns, ['model1_' obj{1}])));
best2 = together.x(find(strcmp(mergedModelTogether.rxns, ['model2_' obj{2}])));

interaction = '';
%determine type of interaction
if best1 < alone1*0.9 %(-)
    if best2 < alone2*0.9 %(-)
        interaction = 'competition';
    elseif best2 > alone2*1.1 %(+)
        interaction = 'parasitism';
    else
        interaction = 'amensalism';
    end
    
elseif best1 > alone1*1.1%(+)
    
    if best2 < alone2*0.9 %(-)
        interaction = 'parasitism';
    elseif best2 > alone2*1.1 %(+)
        interaction = 'mutualism';
    else
        interaction = 'commensalism';
    end
    
else %(N)
    
    if best2 < alone2*0.9 %(-)
        interaction = 'amensalism';
    elseif best2 > alone2*1.1 %(+)
        interaction = 'commensalism';
    else
        interaction = 'neutralism';
    end
end
   

% Competition: (-) and (-)
% Parasitism: one grows faster (+) and the other lower (-)
% Amensalism: one microbe grows slower (-) and the other is unaffected (N)
% Neutralism: (N), (N)
% Commensalism: (+) and (N)
% Mutualism: (+) and (+)

end