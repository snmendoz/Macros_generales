function [interaction, alone1, alone2, best1, best2] = DetermineInteraction(runID, model1, model2, obj, medium, species)

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
[mergedModel, exchangeRxns] = mergeModels(models, shared, 'lbPool', lb, 'indLb', lb_pool, 'indMets', indMets);

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

mergedModelParetoManual = mergedModel;
xs = [];
ys = [];
for i = 0:0.01:1
    alpha = i;
    beta = 1 - alpha;
    mergedModelParetoManual.c(posObj(1)) = alpha;
    mergedModelParetoManual.c(posObj(2)) = beta;
    fba_pm = optimizeCbModel(mergedModelParetoManual);
    xs = [xs; fba_pm.x((posObj(1)))];
    ys = [ys; fba_pm.x((posObj(2)))];
end


%run pareto
[x, y, x_adicionales, y_adicionales] = runBensolve(runID, mergedModel, posObj , 'output', 1);
f=figure('Visible','off');
plot(xs, ys);
hold on
plot(x,y, 'r.')
% plot(x_adicionales, y_adicionales)
title(['Pareto front for growth rates of ' species{1} ' and ' species{2} ],'FontSize', 20)
ylabel([species{2} ' growth rate [1/h]'],'FontSize', 18)
xlabel([species{1} ' growth rate [1/h]'],'FontSize', 18)
set(gca,'fontsize',16)
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 33 15]);
saveas(gcf,['Pareto_' species{1} '_' species{2} '.png']);
set(gcf,'PaperOrientation','landscape')
set(gcf,'PaperPosition', [1 1 28 19])
saveas(gcf,['Pareto_' species{1} '_' species{2}  '.pdf'])
close(f);
%determine best pareto point

total = xs + ys;
[maximo, posMaximo] = max(total);

xMaximo = xs(posMaximo); 
best1 = xMaximo;
yMaximo = ys(posMaximo);
best2 = yMaximo;

interaction = '';
%determine type of interaction
if xMaximo < fba1b.f*0.9 %(-)
    if yMaximo < fba2b.f*0.9 %(-)
        interaction = 'competition';
    elseif yMaximo > fba2b.f*1.1 %(+)
        interaction = 'parasitism';
    else
        interaction = 'amensalism';
    end
    
elseif xMaximo > fba1b.f*1.1%(+)
    
    if yMaximo < fba2b.f*0.9 %(-)
        interaction = 'parasitism';
    elseif yMaximo > fba2b.f*1.1 %(+)
        interaction = 'mutualism';
    else
        interaction = 'commensalism';
    end
    
else %(N)
    
    if yMaximo < fba2b.f*0.9 %(-)
        interaction = 'amensalism';
    elseif yMaximo > fba2b.f*1.1 %(+)
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