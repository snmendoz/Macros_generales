function [r_c_analysis,s_p_analysis] = makeReducedCostAnalysis(model,fba, media, pos_ex, fileName,threshold)
%reduced cost. 
% in a maximization problem, the reduced cost of a variable (reaction) is the 
% amount by which the coeficient for that variable in the objective
% function would have to increase to see a positive value (non-zero) in the
%optimal solution

%in a maximization problem is the change in the value of the objective function
% in response to a change in the right hand side of one constraint (mass balance of one particular metabolite)

if nargin <6
    threshold = 10^-9;
end
if isempty(fba)
    fba = optimizeCbModel(model);
end

pos_sp = find(abs(fba.y)>threshold);
sp = fba.y(pos_sp);
[sp_sorted, ind_sp_sorted] = sort(sp,'descend');
mets_sp = model.mets(pos_sp(ind_sp_sorted));
metsNames_sp = model.metNames(pos_sp(ind_sp_sorted));
s_p_analysis = [mets_sp metsNames_sp num2cell(sp_sorted)];

pos_rc = find(abs(fba.w)>threshold);
rc = fba.w(pos_rc);
src = (rc.*fba.x(pos_rc))/fba.f;

[rc_sorted, ind_sorted] = sort(rc,'descend');
src_sorted = src(ind_sorted);
rxns_rc = model.rxns(pos_rc(ind_sorted));
rxnNames_rc = model.rxnNames(pos_rc(ind_sorted));
rxns_lb_rc = model.lb(pos_rc(ind_sorted));
rxns_ub_rc = model.ub(pos_rc(ind_sorted));

fluxes = fba.x(pos_rc(ind_sorted));

%test
model = changeRxnBounds(changeRxnBounds(model, model.rxns(pos_ex), 0, 'l'), media.reactions, media.lb, 'l');
fba_base = optimizeCbModel(model);

real_effect_lb_down = zeros(length(rxns_rc),1);
real_effect_lb_up = zeros(length(rxns_rc),1);
real_effect_ub_down = zeros(length(rxns_rc),1);
real_effect_ub_up = zeros(length(rxns_rc),1);

for i = 1:length(rxns_rc)
    model_test = changeRxnBounds(changeRxnBounds(model, model.rxns(pos_ex), 0, 'l'), media.reactions, media.lb, 'l');
    model_test = changeRxnBounds(model_test,rxns_rc(i),rxns_lb_rc(i)-1,'l');
    fba_test = optimizeCbModel(model_test);
    real_effect_lb_down(i) = fba_test.f/fba_base.f;
    
    if rxns_lb_rc(i)+1>rxns_ub_rc(i)
        real_effect_lb_up(i) = nan;
    else
        model_test = changeRxnBounds(changeRxnBounds(model, model.rxns(pos_ex), 0, 'l'), media.reactions, media.lb, 'l');
        model_test = changeRxnBounds(model_test,rxns_rc(i),rxns_lb_rc(i)+1,'l');
        fba_test = optimizeCbModel(model_test);
        real_effect_lb_up(i) = fba_test.f/fba_base.f;
    end
    
    if rxns_ub_rc(i)-1<rxns_lb_rc(i)
        real_effect_ub_down(i) = nan;
    else
        model_test2 = changeRxnBounds(changeRxnBounds(model, model.rxns(pos_ex), 0, 'l'), media.reactions, media.lb, 'l');
        model_test2 = changeRxnBounds(model_test2,rxns_rc(i),rxns_ub_rc(i)-1,'u');
        fba_test2 = optimizeCbModel(model_test2);
        real_effect_ub_down(i) = fba_test2.f/fba_base.f;
    end
    
    model_test2 = changeRxnBounds(changeRxnBounds(model, model.rxns(pos_ex), 0, 'l'), media.reactions, media.lb, 'l');
    model_test2 = changeRxnBounds(model_test2,rxns_rc(i),rxns_ub_rc(i)+1,'u');
    fba_test2 = optimizeCbModel(model_test2);
    real_effect_ub_up(i) = fba_test2.f/fba_base.f;
    
end

newFluxes = zeros(length(rxns_rc),1);
for i = 1:length(rxns_rc)
    modelAux = model;
    if rc_sorted(i)>0
        new_coef = rc_sorted(i)*1.1;
    else
        new_coef = rc_sorted(i)*1.1;
    end
     
    modelAux.c(getPosOfElementsInArray(rxns_rc(i), model.rxns)) = new_coef;
    fba_test3 = optimizeCbModel(modelAux);
    newFluxes(i) = fba_test3.x(getPosOfElementsInArray(rxns_rc(i), model.rxns)); 
end

labels = {'rxns','rxnNames','rxnEquations','rxnEquations','lowerBound','upperBound','flux','newFluxes','lowerBound_down',...
    'lowerBound_up','upperBound_down','lowerBound_up','reducedCost','ScaledReducedCost'};
r_c_analysis = [rxns_rc, rxnNames_rc, getRxn_cobraFormat(model, rxns_rc),getRxn_cobraFormat(model, rxns_rc,1),num2cell(rxns_lb_rc),...
    num2cell(rxns_ub_rc), num2cell(fluxes), num2cell(newFluxes),num2cell(real_effect_lb_down), num2cell(real_effect_lb_up), num2cell(real_effect_ub_down), num2cell(real_effect_ub_up), num2cell(rc_sorted), num2cell(src_sorted)];
r_c_analysis = [labels;r_c_analysis];
xlswrite(fileName,r_c_analysis)

end