function [required, alternatives, all_sets] = cEMAF(communityModel, constraints, exchangeRxns, subjectToMinFlux, notAllowed)

if nargin < 4
    subjectToMinFlux = 0;
end
if nargin < 5
    notAllowed = {};
end

model = communityModel;

%% first transform your model to irreversible if needed
posExchange = getPosOfElementsInArray(exchangeRxns,model.rxns);

if subjectToMinFlux  
    
    for i = 1:length(constraints.rxnList)
        if strcmp(constraints.sense(i),'G')
            model = changeRxnBounds(model, constraints.rxnList(i), constraints.values(i),'l');
        elseif strcmp(constraints.sense(i),'L')
            model = changeRxnBounds(model, constraints.rxnList(i), constraints.values(i),'u');
        else
            model = changeRxnBounds(model, constraints.rxnList(i), constraints.values(i),'b');
        end
    end
    tol = 10^-6;
    [minimizedFlux, modelIrrev, ~, ~, rev2irrev, ~]= minimizeModelFlux_local(model);
    modelIrrev = changeRxnBounds(modelIrrev,'netFlux',minimizedFlux.f+tol,'u');
    doubleDirection = cellfun(@length, rev2irrev(posExchange))==2;
    posExchangeRxns = posExchange;
    posExchangeRxns(find(doubleDirection==0)) = cell2mat(rev2irrev(posExchange(find(doubleDirection==0))));
    posExchangeRxns(find(doubleDirection==1)) = getPosOfElementsInArray(strcat(exchangeRxns(find(doubleDirection==1)),'_b'), modelIrrev.rxns);
    exchangeRxnsIrrev = modelIrrev.rxns(posExchangeRxns);
else
    if any(model.lb(posExchange)<0)
        [modelIrrev,~,rev2irrev,~] = convertToIrreversible(model)%,'sRxns',exchangeRxns,'OrderReactions',true);
        doubleDirection = cellfun(@length, rev2irrev(posExchange))==2;
        posExchangeRxns = posExchange;
        posExchangeRxns(find(doubleDirection==0)) = cell2mat(rev2irrev(posExchange(find(doubleDirection==0))));
        posExchangeRxns(find(doubleDirection==1)) = getPosOfElementsInArray(strcat(exchangeRxns(find(doubleDirection==1)),'_b'), modelIrrev.rxns);
        exchangeRxnsIrrev = modelIrrev.rxns(posExchangeRxns);
    else
        modelIrrev = model;
        posExchangeRxns = posExchange;
        exchangeRxnsIrrev = modelIrrev.rxns(posExchangeRxns);
    end
end

%% Write MILP problem
H = 1000;

S = modelIrrev.S;
[m,n_rxns] = size(S);
n_rxns_int = length(posExchangeRxns);
n_variables_total = n_rxns + n_rxns_int;


% Helper arrays for extracting solutions
sel_cont_sol = 1:n_rxns;
sel_int_sol = n_rxns+1:n_rxns+n_rxns_int;

% Set variable types

vartype_bl(1:n_variables_total) = 'C';
vartype_bl(n_rxns+1:n_rxns+n_rxns_int) = 'B';

% Set bilevel objective
c_bl = zeros(n_variables_total,1);
c_bl(n_rxns+1:n_rxns+n_rxns_int,1) = ones(n_rxns_int,1);

% S*v = 0
A_bl = [S zeros(m,n_rxns_int)];
b_bl = modelIrrev.b;
csense_bl = modelIrrev.csense;

%lb
lb_bl = zeros(n_variables_total,1);
lb_bl(1:n_rxns) = modelIrrev.lb;
lb_bl(vartype_bl == 'B') = 0;

%ub
ub_bl = zeros(n_variables_total,1);
ub_bl(1:n_rxns) = modelIrrev.ub;
ub_bl(vartype_bl == 'B') = 1;

%constraints

pos_rxn_constraints = getPosOfElementsInArray(constraints.rxnList,modelIrrev.rxns);

A_constraints = zeros(length(constraints.rxnList),n_variables_total);
b_constraints = zeros(length(constraints.rxnList),1);
csense_constraints = repmat('E',length(constraints.rxnList),1);

for i = 1:length(pos_rxn_constraints)
    A_constraints(i,pos_rxn_constraints(i)) = 1;
    b_constraints(i) = constraints.values(i);
    csense_constraints(i) = constraints.sense(i);
end

A_bl = [A_bl;A_constraints];
b_bl = [b_bl;b_constraints];
csense_bl = [csense_bl;csense_constraints];

% v(j) <= ub(j)*y(j) j in selected
% v(j) - ub(j)*y(j) <= 0 in selected
A_int_constraints = zeros(n_rxns_int,n_variables_total);
for i = 1:n_rxns_int
    A_int_constraints(i,posExchangeRxns(i)) = H;
    A_int_constraints(i,n_rxns+i) = -H*modelIrrev.ub(posExchangeRxns(i));
end

A_bl = [A_bl; A_int_constraints];
b_bl = [b_bl; zeros(n_rxns_int,1)];
csense_bl(end+1:end+n_rxns_int) = 'L';

% combinations not allowed

n_notAllowed = 0;
for i = 1:size(notAllowed,1)
    for j = 1:size(notAllowed,2)
        if i>=j; continue; end
        n_notAllowed = n_notAllowed+length(notAllowed{i,j});
    end
end
if n_notAllowed>0
    A_notAllowed = zeros(n_notAllowed,n_variables_total);
    b_notAllowed = zeros(n_notAllowed,1);
    csense_notAllowed = repmat('E',n_notAllowed,1);
    cont = 0;
    for i = 1:size(notAllowed,1)
        for j = 1:size(notAllowed,2)
            if i>=j; continue; end
            
            %species j will not provide to species i the following metabolites
            %that means that the system exchange should be at least as low as
            %the uptake of species i.
            %v_system_exchange_met_k  >= uptake_met_k_species_i
            %v_system_exchange_met_k - uptake_met_k_species_i >= 0
            n_notAllowed = n_notAllowed+length(notAllowed{i,j});
            for k = 1:length(notAllowed{i,j})
                cont = cont+1;
                pos_sys_exc_met_k = getPosOfElementsInArray({['EX_',notAllowed{i,j}{k},'_pool_b']}, modelIrrev.rxns);
                pos_uptake_met_k = getPosOfElementsInArray({['model',num2str(i),'_EX_',notAllowed{i,j}{k},'_e_b']}, modelIrrev.rxns);
                A_notAllowed(cont,pos_sys_exc_met_k) = 1;
                A_notAllowed(cont,pos_uptake_met_k) = -1;
                csense_notAllowed(cont) = 'G';
            end
        end
    end
    A_bl = [A_bl; A_notAllowed];
    b_bl = [b_bl; b_notAllowed];
    csense_bl = [csense_bl;csense_notAllowed];
end

% Construct problem structure
bilevelMILPproblem.A = A_bl;
bilevelMILPproblem.b = b_bl;
bilevelMILPproblem.c = c_bl;
bilevelMILPproblem.csense = csense_bl;
bilevelMILPproblem.lb = lb_bl;
bilevelMILPproblem.ub = ub_bl;
bilevelMILPproblem.vartype = vartype_bl;
bilevelMILPproblem.contSolInd = sel_cont_sol;
bilevelMILPproblem.intSolInd = sel_int_sol;
bilevelMILPproblem.osense = 1;

%% solve MILP
tic
solution = solveCobraMILP(bilevelMILPproblem);
toc

%% enumerate

solutions = {};
solutions_rxns = {};
solutions{1} = findSelectedInSolution(solution);
solutions_rxns{1} = exchangeRxnsIrrev(solutions{1});
lengthSolutions = length(solutions{1});

A_bl = [A_bl; [zeros(1,n_rxns), ones(1,n_rxns_int)]];
b_bl = [b_bl;lengthSolutions];
csense_bl = [csense_bl;'E'];

while strcmp(solution.origStat,'OPTIMAL') && length(findSelectedInSolution(solution))==lengthSolutions
    pos_bin_variables = sel_int_sol(solutions{end});
    
    A_newSolution = zeros(1,n_rxns+n_rxns_int);
    A_newSolution(pos_bin_variables) = ones(1,length(pos_bin_variables));
    
    A_bl = [A_bl; A_newSolution];
    b_bl = [b_bl;lengthSolutions-1];
    csense_bl = [csense_bl;'L'];
    
    bilevelMILPproblem.A = A_bl;
    bilevelMILPproblem.b = b_bl;
    bilevelMILPproblem.csense = csense_bl;
    solution = solveCobraMILP(bilevelMILPproblem);
    if strcmp(solution.origStat,'OPTIMAL') && length(findSelectedInSolution(solution))==lengthSolutions
        solutions{end+1} = findSelectedInSolution(solution);
        save('solutions_EMAF','solutions')
        
        solutions_rxns{end+1} = exchangeRxnsIrrev(solutions{end});
        save('solutions_rxns_EMAF','solutions_rxns')
        
        disp('new solution found')
        disp(length(solutions))
    end
    
end

%% find required and alternative
required = solutions_rxns{1};
for i = 2:length(solutions_rxns)
    required = intersect(required, solutions_rxns{i});
end

rests = cell(length(solutions_rxns),1);
for i = 1:length(solutions_rxns)
    rests{i} = setdiff(solutions_rxns{i}, required);
end

n_alternative_groups = length(solutions_rxns{1})-length(required);
alternatives = cell(n_alternative_groups,1);

n = 1;
for i = 1:length(rests{n})
    rest_i = rests{n}(i);
    others = rests{n}(setdiff(1:length(rests{n}),i));
    alternatives_i = cellfun(@(y) setdiff(y,others), rests(find(cellfun(@(x) all(ismember(others, x)), rests))));
    alternatives{i} = union(rest_i, alternatives_i);
end
alternatives = cellfun(@(x) strjoin(x,', '), alternatives, 'UniformOutput', false);

all_sets = solutions_rxns;

end


function [minimizedFlux, modelIrrev, solution_rev, matchRev, rev2irrev, irrev2rev] = minimizeModelFlux_local(model,GeneOption)
% This function finds the minimum flux through the network and returns the
% minimized flux and an irreversible model convert model to irrev


[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);
% add pseudo-metabolite to measure flux through network

if nargin==1
    GeneOption=0;
end

%Add the metabolite
modelIrrev = addMetabolite(modelIrrev,'fluxMeasure');
%Then set all the Stoichiometric entries.

if GeneOption==0 % signal that you want to minimize the sum of all gene and non-gene associated fluxes
    modelIrrev.S(end,:) = ones(size(modelIrrev.S(1,:)));
elseif GeneOption==1 % signal that you want to minimize the sum of only gene-associated fluxes
    %find all reactions which are gene associated
    Ind=find(sum(modelIrrev.rxnGeneMat,2)>0);
    modelIrrev.S(end,:) = zeros(size(modelIrrev.S(1,:)));
    modelIrrev.S(end,Ind) = 1;
elseif GeneOption==2 % signal that you want to minimize the sum of only NON gene-associated fluxes
    %find all reactions which are gene associated
    Ind=find(sum(modelIrrev.rxnGeneMat,2)==0);
    modelIrrev.S(end,:) = zeros(size(modelIrrev.S(1,:)));
    modelIrrev.S(end,Ind) = 1;
end

% add a pseudo reaction that measures the flux through the network
%     modelIrrev = addReaction(modelIrrev,'netFlux',{'fluxMeasure'},[-1],false,0,inf,0,'','');
S_new_rxn = zeros(size(modelIrrev.S,1),1);
S_new_rxn(getPosOfElementsInArray({'fluxMeasure'},modelIrrev.mets),1) = -1;
modelIrrev.S = [modelIrrev.S, sparse(S_new_rxn)];
clear S_new_rxn
modelIrrev.rxns = [modelIrrev.rxns;{'netFlux'}];
modelIrrev.rxnNames = [modelIrrev.rxnNames;{'netFlux'}];
modelIrrev.subSystems = [modelIrrev.subSystems; {''}];
modelIrrev.lb = [modelIrrev.lb; 0];
modelIrrev.ub = [modelIrrev.ub; 1000];
modelIrrev.grRules = [modelIrrev.grRules; {''}];
%     modelIrrev.rules = [modelIrrev.rules; {''}];
%     modelIrrev.rxnGeneMat = [modelIrrev.rxnGeneMat;zeros(1,length(modelIrrev.genes))];
modelIrrev.match = [modelIrrev.match; 0];

% set the flux measuring demand as the objective
modelIrrev.c = zeros(length(modelIrrev.rxns),1);
modelIrrev = changeObjective(modelIrrev, 'netFlux');

% minimize the flux measuring demand (netFlux)
minimizedFlux = optimizeCbModel(modelIrrev,'min');
solution_irrev = minimizedFlux.x;
solution_rev = zeros(size(model.rxns));
solution_rev(irrev2rev(find(matchRev==0))) = solution_irrev(find(matchRev==0));
solution_rev(find(cellfun(@length, rev2irrev)==2)) = cellfun(@(x) solution_irrev(x(1))-solution_irrev(x(2)) ,rev2irrev(find(cellfun(@length, rev2irrev)==2)));

pos_irrev = find(cellfun(@length, rev2irrev)==1);
pos_different = pos_irrev(find(strcmp(model.rxns(pos_irrev), modelIrrev.rxns(pos_irrev)) ==0));
solution_rev(pos_different) = -solution_rev(pos_different);

end

function pos_selected = findSelectedInSolution(solution)

tol = 10^-1;
pos_selected = find(solution.int>1-tol);

end