function [limitingNutrients, improvement] = findLimitingNutrientBruteForce(model)

[pos_excRxns, excRxns] = findExcRxnsWithIDs(model, 1);
improve_when_more_available = zeros(size(excRxns));
improvement = zeros(size(excRxns));
increase_parameter = 0.01;

fba = optimizeCbModel(model);

for i = 1:length(excRxns)
    fbaAux = optimizeCbModel(changeRxnBounds(model,excRxns{i} ,model.lb(pos_excRxns(i))-increase_parameter,'l'));
    if fbaAux.f>fba.f+10^-8 || fbaAux.f>fba.f*1.001
        improvement(i) = fbaAux.f/fba.f;
        improve_when_more_available(i) = 1;
    end
end

limitingNutrients = excRxns(improve_when_more_available==1);
improvement = improvement(improve_when_more_available==1);

end