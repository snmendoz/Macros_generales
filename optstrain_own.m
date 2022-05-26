function [modelOutputOptStrain, solutions_rxns, solutions, lengthSolutions, translatedSolutions, timePerSolution, runTimes, n_binary] = optstrain_own(model, product, database, constraints)

%check that database and models are complaint (compartment names)
if allMetsInCOBRAFormat(database.mets) && allMetsInCBMPYFormat(model.mets)
    database = transformModelToCBMPYFormat(database);
end
if allMetsInCBMPYFormat(database.mets) && allMetsInCOBRAFormat(model.mets)
    database = transformModelToCOBRAFormat(database);
end

%get candidates
candidates = database.rxns;
eqs_candidates = getRxn_cobraFormat(database, candidates);

[rxnFormulaWasFound, posRxn, ~, matchWithoutRemovingProtons, ...
    ~, posRxnsFullMatch, posRxnsPartialMatch, validPartialMatchs] = ...
    cellfun(@(x) reactionFormulaInModel(model, x, 1, 1), eqs_candidates,'UniformOutput',0);

pos_realCandidates = ~cell2mat(rxnFormulaWasFound);
candidates = candidates(pos_realCandidates);
eqs_candidates = eqs_candidates(pos_realCandidates);
modelInputOptStrain = model;
%add Reactions From Database
for i = 1:length(candidates)
    modelInputOptStrain = addReaction(modelInputOptStrain,candidates{i},'reactionFormula',eqs_candidates{i},'printLevel',0);
end
%convert reactions to irreversible
[modelInputOptStrain_irrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(modelInputOptStrain,'sRxns',candidates,'OrderReactions',true);

candidatesToRemove_irrev_pos = [];
for i = 1:length(candidates)
    pos_i = getPosOfElementsInArray(candidates(i),modelInputOptStrain.rxns);
    candidatesToRemove_irrev_pos = union(candidatesToRemove_irrev_pos,rev2irrev{pos_i});
end
candidates_irrev = modelInputOptStrain_irrev.rxns(candidatesToRemove_irrev_pos);

n_binary = length(candidates_irrev);

bilevelMILPproblem = buildMILPproblem_gapfilling_sparse({modelInputOptStrain_irrev}, candidates_irrev, {constraints}, [],{},[]);
tic
solution = solveCobraMILP(bilevelMILPproblem);
runTimes = toc;

if strcmp(solution.origStat,'INFEASIBLE')
    solutions = {};
    solutions_rxns = {};
    translatedSolutions =solutions_rxns;
    modelOutputOptStrain = [];
    lengthSolutions= 0;
    timePerSolution = mean(runTimes);
else
    
    solutions = {};
    solutions_rxns = {};
    solutions{1} = findSelectedInSolution(solution);
    solutions_rxns{1} = bilevelMILPproblem.int_ids(solutions{1});
    lengthSolutions = length(solutions{1});
    minValue = solution.obj;
    
    notRemoved = candidates_irrev(findSelectedInSolution(solution));
    removed = setdiff(candidates_irrev, notRemoved);
    pos_notZeros = find(solution.full(getPosOfElementsInArray(removed, modelInputOptStrain_irrev.rxns)));
    removed_flux_notZero = removed(pos_notZeros);
    flux = solution.full(getPosOfElementsInArray(removed_flux_notZero, modelInputOptStrain_irrev.rxns));
    
    modelOutputOptStrain = removeRxns(modelInputOptStrain_irrev,removed);
    
    fba_check = optimizeCbModel(modelOutputOptStrain);
    disp('')
    
    %%
    while strcmp(solution.origStat,'OPTIMAL') && length(findSelectedInSolution(solution))<=lengthSolutions
        
        tic
        bilevelMILPproblem = buildMILPproblem_gapfilling_sparse({modelInputOptStrain_irrev}, candidates_irrev, {constraints}, solutions_rxns,{},[]);
        time = toc;
        runTimes = [runTimes; time];
        solution = solveCobraMILP(bilevelMILPproblem);
        
        if strcmp(solution.origStat,'OPTIMAL') && length(findSelectedInSolution(solution))<=lengthSolutions
            solutions{end+1} = findSelectedInSolution(solution);
            solutions_rxns{end+1} = bilevelMILPproblem.int_ids(solutions{end});
            
            lengthSolutions = min(length(solutions{end}), lengthSolutions);
            
            disp('new solution found')
            disp(length(solutions))
        end
    end
    
    % solutions_info = []
    % for k = 2:length(solutions_rxns)
    % solutions_info = [solutions_info; solutions_rxns{k}'];
    % end
    
    timePerSolution = mean(runTimes);
    
    translatedSolutions =solutions_rxns;
    for i = 1:length(translatedSolutions)
        translatedSolutions{i} = regexprep(translatedSolutions{i}, {'_b$','_f$','_r$'},{'','',''});
    end
    
end



end

function pos_selected = findSelectedInSolution(solution)

tol = 10^-1;
pos_selected = find(solution.int>1-tol);

end