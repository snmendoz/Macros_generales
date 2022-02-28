function [improveWhenRemoved, ratios_removed, improveWhenAdded, ratios_added] = costAnalysis(model,filterout)

fbaBasal = optimizeCbModel(model);
improveWhenRemoved = {};
ratios_removed = [];
improveWhenAdded = {};
ratios_added = [];

mets = getMetsInCompartment(model, 'c');
mets = setdiff(mets,filterout);
for i = 1:length(mets)
    if ismember(mets{i},{'dhmpt_c','dhlam_c','dhf_c'})
        disp(i)
    end
    [modelaux1,rxns] = addDemandReaction(model, mets{i});
    modelaux1 = changeRxnBounds(modelaux1, rxns, 0,'l');
    modelaux1 = changeRxnBounds(modelaux1, rxns,0.01,'u');
    fba1 = optimizeCbModel(modelaux1);

    modelaux2 = changeRxnBounds(modelaux1, rxns,-0.01,'l');
    modelaux2 = changeRxnBounds(modelaux2, rxns,0,'u');
    fba2 = optimizeCbModel(modelaux2);
    
    if fba1.f>fbaBasal.f
        improveWhenRemoved = [improveWhenRemoved;mets(i)];
        ratios_removed = [ratios_removed;fba1.f/fbaBasal.f];
    end
    
    if fba2.f>fbaBasal.f
        improveWhenAdded = [improveWhenAdded;mets(i)];
        ratios_added = [ratios_added;fba2.f/fbaBasal.f];
    end
    
end

end