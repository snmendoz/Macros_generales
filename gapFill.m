function [modelGapFilled, solutions_rxns, solutions, lengthSolutions, translatedSolutions] = gapFill(model, database, constraints, penalizeGeneAssociated, penalizationFactor)

modelGapFilled = [];

if nargin<4 || isempty(penalizeGeneAssociated)
    penalizeGeneAssociated = 0;
end
if nargin<5 || isempty(penalizationFactor)
    penalizationFactor = 3;
end


modelOut = model;

%check that 

%check that the model and the database are consistent (common reactions are equal)
shared = intersect(modelOut.rxns, database.rxns);


[rxnFormulaWasFound, posRxn, ~, matchWithoutRemovingProtons, ...
            ~, posRxnsFullMatch, posRxnsPartialMatch, validPartialMatchs] = ...
            cellfun(@(x) reactionFormulaInModel(database, x, 1, 1), getRxn_cobraFormat(modelOut, shared),'UniformOutput',0);

%rename those which are equal in name but different in reaction equation        
pos_different = find((cell2mat(rxnFormulaWasFound))==0);
for i = 1:length(pos_different)
    pos = getPosOfElementsInArray(shared(pos_different(i)), modelOut.rxns);
    modelOut.rxns{pos} = [modelOut.rxns{pos} '_2'];
end

%make sure that by adding all the reactions, you get growth
diff = setdiff(database.rxns,modelOut.rxns);
fba_check_before = optimizeCbModel(modelOut);

% modelaux = modelOut;
% modelaux = addReactionFromModelRef(modelaux, diff, database);
% fba_check_after = optimizeCbModel(modelaux);

candidatesToRemove = diff;
%exclude exchange reactions
candidatesToRemove = setdiff(candidatesToRemove, database.rxns(findExcRxnsWithIDs(database))); 

fba_check_database = optimizeCbModel(database);

%exclude essential reactions
[essentials,gr,pos_essential] = findEssentialRxns(database, fba_check_database.f*0.8, candidatesToRemove);
candidatesToRemove = setdiff(candidatesToRemove,candidatesToRemove(essentials));
%determine blocked reactions. That will be removed

oldTol= 1e-9;
[minimums1, maximums1] = fluxVariability(database,0,'max',candidatesToRemove);
pos_blocked_1 = intersect(find(minimums1>-oldTol), find(maximums1<oldTol));
blockedRxns = candidatesToRemove(pos_blocked_1);
database_wo_blocked_rxns = removeRxns(database,blockedRxns);
fba_wo_blocked_rxns = optimizeCbModel(database_wo_blocked_rxns);
if fba_wo_blocked_rxns.f > 1e-3
    candidatesToRemove = setdiff(candidatesToRemove,blockedRxns);
    database = database_wo_blocked_rxns;
end

[databaseIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(database,'sRxns',candidatesToRemove,'OrderReactions',true);
candidatesToRemove_irrev_pos = [];

doubleDirection = zeros(size(candidatesToRemove));
for i = 1:length(candidatesToRemove)
    pos_i = getPosOfElementsInArray(candidatesToRemove(i),database.rxns);
    candidatesToRemove_irrev_pos = union(candidatesToRemove_irrev_pos,rev2irrev{pos_i});
    doubleDirection(i) = length(rev2irrev{pos_i})==2;
end
candidatesToRemove_irrev = databaseIrrev.rxns(candidatesToRemove_irrev_pos);

% databaseIrrev.rxns = strcat('model1_', databaseIrrev.rxns);
% databaseIrrev.mets = strcat('model1_', databaseIrrev.mets);

if penalizeGeneAssociated
    penalizedReactions = candidatesToRemove_irrev(find(cellfun(@isempty,databaseIrrev.grRules(getPosOfElementsInArray(candidatesToRemove_irrev, databaseIrrev.rxns)))==0));
    penalizationFactors = repmat(penalizationFactor, length(penalizedReactions),1);
else
    penalizedReactions = {};
    penalizationFactors = [];
end
bilevelMILPproblem = buildMILPproblem_gapfilling({databaseIrrev}, candidatesToRemove_irrev, {constraints}, [], penalizedReactions, penalizationFactors);
tic
solution = solveCobraMILP(bilevelMILPproblem);
toc

% {'ADCL'    }
%     {'ADCS'    }
%     {'CHORS'   }
%     {'DDPA'    }
%     {'DHQS'    }
%     {'DHQTi_f' }
%     {'PSCVT_f' }
%     {'SHK3Dr_f'}
%     {'SHKK'    }
%     {'TKT2_b'  }
%% enumerate
solutions = {};
solutions_rxns = {};
solutions{1} = findSelectedInSolution(solution);
solutions_rxns{1} = bilevelMILPproblem.int_ids(solutions{1});
lengthSolutions = length(solutions{1});
minValue = solution.obj;

notRemoved = candidatesToRemove_irrev(findSelectedInSolution(solution));
removed = setdiff(candidatesToRemove_irrev, notRemoved);
pos_notZeros = find(solution.full(getPosOfElementsInArray(removed, databaseIrrev.rxns)));
removed_flux_notZero = removed(pos_notZeros);
flux = solution.full(getPosOfElementsInArray(removed_flux_notZero, databaseIrrev.rxns));

fba_check = optimizeCbModel(removeRxns(databaseIrrev,removed))
disp('')

%%

% getPosOfElementsInArray(removed_irrev(find(solution.cont(getPosOfElementsInArray(removed_irrev, databaseIrrev.rxns)'))), candidatesToRemove_irrev)
% sparse(bilevelMILPproblem.A(length(databaseIrrev.mets)+1+21,:))
% databaseIrrev.ub(46)*1e6
% bilevelMILPproblem.A(length(databaseIrrev.mets)+1+21,46)*solution.full(46)+bilevelMILPproblem.A(length(databaseIrrev.mets)+1+21,922)*solution.full(922)
% 
% databaseIrrev.lb(46)*1e6
% sparse(bilevelMILPproblem.A(length(databaseIrrev.mets)+1+length(candidatesToRemove_irrev)+21,:))
% bilevelMILPproblem.A(length(databaseIrrev.mets)+1+length(candidatesToRemove_irrev)+21,46)*solution.full(46)

% sparse(bilevelMILPproblem.A(length(databaseIrrev.mets)+1+length(candidatesToRemove_irrev)+377,:))
% bilevelMILPproblem.A(length(databaseIrrev.mets)+1+length(candidatesToRemove_irrev)+377,777)*solution.full(777)

% sparse(bilevelMILPproblem.A(length(databaseIrrev.mets)+1+377,:))
% bilevelMILPproblem.A(length(databaseIrrev.mets)+1+377,777)*solution.full(777)+bilevelMILPproblem.A(length(databaseIrrev.mets)+1+377,1278)*solution.full(1278)

% sparse(bilevelMILPproblem.A(length(databaseIrrev.mets)+1+11,:))
% bilevelMILPproblem.A(length(databaseIrrev.mets)+1+11,28)*solution.full(28)+bilevelMILPproblem.A(length(databaseIrrev.mets)+1+11,1278)*solution.full(28)
%%
while strcmp(solution.origStat,'OPTIMAL') && length(findSelectedInSolution(solution))==lengthSolutions
    
    bilevelMILPproblem = buildMILPproblem_gapfilling({databaseIrrev}, candidatesToRemove_irrev, {constraints}, solutions_rxns,penalizedReactions, penalizationFactors);
    
    solution = solveCobraMILP(bilevelMILPproblem);
    
    if strcmp(solution.origStat,'OPTIMAL') && length(findSelectedInSolution(solution))==lengthSolutions
        solutions{end+1} = findSelectedInSolution(solution);
        save('solutions','solutions')
        
        solutions_rxns{end+1} = bilevelMILPproblem.int_ids(solutions{end});
        save('solutions_rxns','solutions_rxns')
        
        disp('new solution found')
        disp(length(solutions))
    end
end


% translate solutions
translatedSolutions = cell(size(solutions_rxns));
removedRxns = cell(size(solutions_rxns));
for i = 1:length(translatedSolutions)
    selected_rxns_irrev = regexprep(solutions_rxns{i}, 'bin_','');
    removed_irrev = setdiff(candidatesToRemove_irrev,selected_rxns_irrev);
    
    pos_removed = cellfun(@(x) all(ismember(databaseIrrev.rxns(rev2irrev{getPosOfElementsInArray({x},database.rxns)}), removed_irrev)), candidatesToRemove);
    removed = candidatesToRemove(pos_removed);
    removedRxns{i} = removed;
    notRemoved = setdiff(candidatesToRemove,removed);
    translatedSolutions{i} = notRemoved;
    
%     keepAsItWasOriginally = 
%     keepOnlyBackward = 
%     keepOnlyForward =
end
    
if length(solutions)>=1
    modelGapFilled = removeRxns(database, removedRxns{1});
    %remove unneeded exchangeReactions
    [metConn, ~] = networkTopology(modelGapFilled);
    orphansRxns_wrtMets = find(arrayfun(@(x) all(ismember(find(modelGapFilled.S(:, x)), find(metConn==1))), 1:length(modelGapFilled.rxns)));
    modelGapFilled = removeRxns(modelGapFilled, intersect(modelGapFilled.rxns(orphansRxns_wrtMets), diff)); 
end
end

function pos_selected = findSelectedInSolution(solution)

tol = 10^-1;
pos_selected = find(solution.int>1-tol);

end