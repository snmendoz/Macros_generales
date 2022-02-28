function [minimal_set, impact, order, impact_order] = makeSensitivityAnalysis_DemandReactionsApproach(model, direction ,leftOutMets, approach)

factor_forward = 1.01;
factor_backward = 0.99;
order = {};
impact_order = {};

switch approach
    
    case 'forward'
        minimal_set = {};
        specific_mets = setdiff(model.mets, leftOutMets);
        fba = optimizeCbModel(model);
        for i = 1:length(specific_mets)
            [model_aux, dm_rxn] = addDemandReaction(model, specific_mets(i));
            if strcmp(direction, 'forward')
                model_aux = changeRxnBounds(model_aux,dm_rxn,0,'l');
                model_aux = changeRxnBounds(model_aux,dm_rxn,1000,'u');
            elseif strcmp(direction, 'backward')
                model_aux = changeRxnBounds(model_aux,dm_rxn,-1000,'l');
                model_aux = changeRxnBounds(model_aux,dm_rxn,0,'u');
            elseif strcmp(direction, 'both')
                model_aux = changeRxnBounds(model_aux,dm_rxn,-1000,'l');
                model_aux = changeRxnBounds(model_aux,dm_rxn,1000,'u');
            end
            fba_aux = optimizeCbModel(model_aux);
            if fba_aux.f > fba.f*factor_forward
                minimal_set = union(minimal_set, dm_rxn);
            end
        end
        impact = zeros(length(minimal_set),1);
        for i = 1:length(minimal_set)
            [model_aux, dm_rxn] = addDemandReaction(model, regexprep(minimal_set{i}, '^DM_', ''));
            if strcmp(direction, 'forward')
                model_aux = changeRxnBounds(model_aux,dm_rxn,0,'l');
                model_aux = changeRxnBounds(model_aux,dm_rxn,1000,'u');
            elseif strcmp(direction, 'backward')
                model_aux = changeRxnBounds(model_aux,dm_rxn,-1000,'l');
                model_aux = changeRxnBounds(model_aux,dm_rxn,0,'u');
            elseif strcmp(direction, 'both')
                model_aux = changeRxnBounds(model_aux,dm_rxn,-1000,'l');
                model_aux = changeRxnBounds(model_aux,dm_rxn,1000,'u');
            end
            model_aux = changeRxnBounds(model_aux, dm_rxn,0.1,'u');
            fba_aux = optimizeCbModel(model_aux);
            impact(i) = fba_aux.f/fba.f;
        end
        
    case 'backward'
        fba_base = optimizeCbModel(model);
        specific_mets = setdiff(model.mets, leftOutMets);
        model_aux = model;
        model_aux = addDemandReaction(model_aux, specific_mets);
        
        demand_rxns = setdiff(model_aux.rxns,model.rxns);
        for i = 1:length(demand_rxns)
            if strcmp(direction, 'forward')
                model_aux = changeRxnBounds(model_aux,demand_rxns{i},0,'l');
                model_aux = changeRxnBounds(model_aux,demand_rxns{i},1000,'u');
            elseif strcmp(direction, 'backward')
                model_aux = changeRxnBounds(model_aux,demand_rxns{i},-1000,'l');
                model_aux = changeRxnBounds(model_aux,demand_rxns{i},0,'u');
            elseif strcmp(direction, 'both')
                model_aux = changeRxnBounds(model_aux,demand_rxns{i},-1000,'l');
                model_aux = changeRxnBounds(model_aux,demand_rxns{i},1000,'u');
            end
        end
        fba = optimizeCbModel(model_aux);
        minimal_set={};
        for i = 1:length(demand_rxns)
            model_aux2 = changeRxnBounds(model_aux,demand_rxns(i),0,'b');
            fba_aux = optimizeCbModel(model_aux2);
            if fba_aux.f<fba.f*factor_backward
                minimal_set = union(minimal_set, demand_rxns(i));
            end
            
            impact = zeros(length(minimal_set),1);
            for i = 1:length(minimal_set)
                model_aux2 = changeRxnBounds(model_aux,minimal_set(i),0,'b');
                fba_aux = optimizeCbModel(model_aux2);
                impact(i) = fba_aux.f/fba.f;
            end
        end
    case 'backward-forward'
        [minimal_set, impact] = makeSensitivityAnalysis_DemandReactionsApproach(model, direction ,leftOutMets, 'backward');
        if length(minimal_set)>0
            order = {};
            impact_order = {};
            all_dm = {};
            excluded =leftOutMets;
            cont = 0;
            model_recursive = model;
            while length(minimal_set)>length(all_dm)
                [minimal_set_aux, impact_aux] = makeSensitivityAnalysis_DemandReactionsApproach(model_recursive, direction, excluded, 'forward');
                if isempty(minimal_set_aux)
                    break;
                end
                [model_recursive, dm] = addDemandReaction(model_recursive, regexprep(minimal_set_aux, '^DM_',''));
                if strcmp(direction, 'backward')
                    model_recursive = changeRxnBounds(model_recursive, dm,-1000,'l');
                    model_recursive = changeRxnBounds(model_recursive, dm,0,'u');
                end
                all_dm = union(all_dm,minimal_set_aux);
                excluded = union(excluded,regexprep(minimal_set_aux, '^DM_',''));
                cont = cont+1;
                order{cont} = minimal_set_aux;
                impact_order{cont} = impact_aux;
            end
            
        end
end
end

