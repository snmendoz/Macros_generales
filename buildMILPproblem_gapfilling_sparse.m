function bilevelMILPproblem = buildMILPproblem_gapfilling_sparse(models, binary_variable_ids, constraints, previousSolutions, penalizedReactions, penalizationFactors)

if nargin < 4
    previousSolutions = [];
end

H1 = 1e-3;
H2 = 1e9;


cont_ids = [];
metNames = [];
for i = 1:length(models)    
    n_rxns_per_model(i) = length(models{i}.rxns);
    n_mets_per_model(i) = length(models{i}.mets);
    cont_ids = [cont_ids;strcat('model',num2str(i),'_',models{i}.rxns)];
    metNames = [metNames;strcat('model',num2str(i),'_',models{i}.mets)];
end
int_ids = strcat('bin_',binary_variable_ids);
variable_ids = [cont_ids;int_ids];

n_rxn_total = sum(n_rxns_per_model);
n_met_total = sum(n_mets_per_model);

n_rxns_int = length(binary_variable_ids);
n_variables_total = n_rxn_total + n_rxns_int;


% Helper arrays for extracting solutions
sel_cont_sol = [1:n_rxn_total];
sel_int_sol = [n_rxn_total+1:n_rxn_total+n_rxns_int];

% Set variable types

vartype_bl(1:n_variables_total) = 'C';
vartype_bl(n_rxn_total+1:n_rxn_total+n_rxns_int) = 'B';

% Set bilevel objective
c_bl = zeros(n_variables_total,1);
c_bl(n_rxn_total+1:n_rxn_total+n_rxns_int,1) = ones(n_rxns_int,1);

if ~isempty(penalizedReactions)
    posPenalizedReactions = getPosOfElementsInArray(penalizedReactions, binary_variable_ids);
    c_bl(n_rxn_total+posPenalizedReactions,1) = penalizationFactors';
end


% S*v = 0
A_bl = sparse(models{1}.S);
for i = 2:length(models)
    A_bl = [[A_bl, sparse(zeros(size(A_bl,1),size(models{i}.S,2))) ]; [sparse(zeros(size(models{i}.S,1),size(A_bl,2))) models{i}.S]];
end
A_bl = [A_bl zeros(size(A_bl,1),n_rxns_int)];
b_bl = zeros(n_met_total,1);
csense_bl(1:n_met_total) = 'E';

%lb
lb = [];
for i = 1:length(models)
    lb = [lb;models{i}.lb];
end
lb_bl(1:n_rxn_total) = lb;
lb_bl(vartype_bl == 'B') = 0;

%ub
ub_bl = zeros(n_variables_total,1);
ub = [];
for i = 1:length(models)
    ub = [ub;models{i}.ub];
end
ub_bl(1:n_rxn_total) = ub;
ub_bl(vartype_bl == 'B') = 1;

%constraints
for k = 1:length(models)  
    constrOpt = constraints{k};
    pos_rxn_constraints = getPosOfElementsInArray(strcat('model',num2str(k),'_',constrOpt.rxnList),variable_ids);
    
    A_constraints = zeros(length(constrOpt.rxnList),n_variables_total);
    b_constraints = zeros(length(constrOpt.rxnList),1);
    csense_constraints = repmat('E',length(constrOpt.rxnList),1);
    
    for i = 1:length(pos_rxn_constraints)
        A_constraints(i,pos_rxn_constraints(i)) = 1;
        b_constraints(i) = constrOpt.values(i);
        csense_constraints(i) = constrOpt.sense(i);
    end
    
    if size(csense_constraints,1)>size(csense_constraints,2)
        csense_constraints = csense_constraints';
    end
    
    A_bl = [A_bl;A_constraints];
    b_bl = [b_bl;b_constraints];
    csense_bl = [csense_bl,csense_constraints];
    
end

% v(j) <= ub(j)*y(j) j in selected
% v(j) - ub(j)*y(j) <= 0 in selected
A_int_constraints = sparse(n_rxns_int*length(models),n_variables_total);
cont = 0;
for i = 1:n_rxns_int
    for j = 1:length(models)
        cont = cont+1;
        pos = getPosOfElementsInArray({strcat('model',num2str(j),'_',binary_variable_ids{i})}, variable_ids);
        A_int_constraints(cont,pos) = H1;
        A_int_constraints(cont,n_rxn_total+i) = -H1*ub(pos);
    end
end
A_bl = [A_bl; A_int_constraints];
b_bl = [b_bl; zeros(n_rxns_int*length(models),1)];
csense_bl(end+1:end+n_rxns_int*length(models)) = 'L';

A_int_constraints = sparse(n_rxns_int*length(models),n_variables_total);
cont = 0;
for i = 1:n_rxns_int
    for j = 1:length(models)
        cont = cont+1;
        pos = getPosOfElementsInArray({strcat('model',num2str(j),'_',binary_variable_ids{i})}, variable_ids);
        A_int_constraints(cont,pos) = H2;
        A_int_constraints(cont,n_rxn_total+i) = -H2*ub(pos);
    end
end
A_bl = [A_bl; A_int_constraints];
b_bl = [b_bl; zeros(n_rxns_int*length(models),1)];
csense_bl(end+1:end+n_rxns_int*length(models)) = 'L';

A_int_constraints = sparse(n_rxns_int*length(models),n_variables_total);
cont = 0;
for i = 1:n_rxns_int
    for j = 1:length(models)
        cont = cont+1;
        pos = getPosOfElementsInArray({strcat('model',num2str(j),'_',binary_variable_ids{i})}, variable_ids);
        A_int_constraints(cont,pos) = 1;
        A_int_constraints(cont,n_rxn_total+i) = -ub(pos);
    end
end
A_bl = [A_bl; A_int_constraints];
b_bl = [b_bl; zeros(n_rxns_int*length(models),1)];
csense_bl(end+1:end+n_rxns_int*length(models)) = 'L';


% v(j) >= lb(j)*y(j) j in selected
% v(j) - lb(j)*y(j) >= 0 in selected
A_int_constraints = zeros(n_rxns_int*length(models),n_variables_total);
cont = 0;
for i = 1:n_rxns_int
    for j = 1:length(models)
        cont = cont+1;
        pos = getPosOfElementsInArray({strcat('model',num2str(j),'_',binary_variable_ids{i})}, variable_ids);
        A_int_constraints(cont,pos) = H1;
        A_int_constraints(cont,n_rxn_total+i) = -H1*lb(pos);
    end
end
A_bl = [A_bl; A_int_constraints];
b_bl = [b_bl; zeros(n_rxns_int*length(models),1)];
csense_bl(end+1:end+n_rxns_int*length(models)) = 'G';
% 
A_int_constraints = sparse(n_rxns_int*length(models),n_variables_total);
cont = 0;
for i = 1:n_rxns_int
    for j = 1:length(models)
        cont = cont+1;
        pos = getPosOfElementsInArray({strcat('model',num2str(j),'_',binary_variable_ids{i})}, variable_ids);
        A_int_constraints(cont,pos) = H2;
        A_int_constraints(cont,n_rxn_total+i) = -H2*lb(pos);
    end
end
A_bl = [A_bl; A_int_constraints];
b_bl = [b_bl; zeros(n_rxns_int*length(models),1)];
csense_bl(end+1:end+n_rxns_int*length(models)) = 'G';

A_int_constraints = sparse(n_rxns_int*length(models),n_variables_total);
cont = 0;
for i = 1:n_rxns_int
    for j = 1:length(models)
        cont = cont+1;
        pos = getPosOfElementsInArray({strcat('model',num2str(j),'_',binary_variable_ids{i})}, variable_ids);
        A_int_constraints(cont,pos) = 1;
        A_int_constraints(cont,n_rxn_total+i) = -lb(pos);
    end
end
A_bl = [A_bl; A_int_constraints];
b_bl = [b_bl; zeros(n_rxns_int*length(models),1)];
csense_bl(end+1:end+n_rxns_int*length(models)) = 'G';


if ~isempty(previousSolutions)
      
    for i = 1:length(previousSolutions)
        
        lengthSolutions = length(previousSolutions{i});
        
        solution_i =  previousSolutions{i};
        pos_bin_variables = getPosOfElementsInArray(solution_i,variable_ids);
        
        A_newSolution = sparse(1,n_rxn_total+n_rxns_int);
        A_newSolution(pos_bin_variables) = ones(1,length(pos_bin_variables));
        
        A_bl = [A_bl; A_newSolution];
        b_bl = [b_bl;lengthSolutions-1];
        csense_bl = [csense_bl,'L'];
    end
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
bilevelMILPproblem.variable_ids = variable_ids;
bilevelMILPproblem.cont_ids = cont_ids;
bilevelMILPproblem.int_ids = int_ids;


end