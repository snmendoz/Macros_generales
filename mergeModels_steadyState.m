function communityModel = mergeModels_steadyState(models, medium, sharedMetaboliteIDs, mu)

n_microorganisms = length(models);
for i = 1:n_microorganisms
    biomassRxns{i} = models{i}.rxns{find(models{i}.c)};
end

indMets = {{}, {}};
lb_pool = {[]; []};
for i = 1:length(models)
    if allMetsInCOBRAFormat(models{i}.mets)
        models{i} = transformModelToCBMPYFormat(models{i});
    end
end
communityModel = contatenate_species(models,medium);


description = 'community model';

pos = find(find(communityModel.ub==1000));
communityModel.ub(pos) = 99999*ones(length(pos),1);

pos = find(find(communityModel.lb==-1000));
communityModel.lb(pos) = -99999*ones(length(pos),1);

for i = 1:length(models)
    pos = find(models{i}.ub>=1000);
    models{i}.ub(pos) = 99999*ones(length(pos),1);
    
    pos = find(models{i}.lb<=-1000);
    models{i}.lb(pos) = -99999*ones(length(pos),1);
end

% system metabolites
m_sys = communityModel.mets;
mNames_sys = communityModel.metNames;

% system reactions
r_sys = communityModel.rxns;
rNames_sys = communityModel.rxnNames;
ub_sys = communityModel.ub;
lb_sys = communityModel.lb;
obj_sys = communityModel.c;
subS_sys = communityModel.subSystems;
grRules_sys = communityModel.grRules;

%add biomass variables
IDs_new_variables = {};
for i = 1:n_microorganisms
    IDs_new_variables = [IDs_new_variables; ['biomass_m' num2str(i)]];
end

names_new_variables = IDs_new_variables;
ub_new_variables = 1000*ones(length(IDs_new_variables),1);
lb_new_variables = zeros(length(IDs_new_variables),1);
obj_new_variables = zeros(length(IDs_new_variables),1);
subS_new_variables = repmat({'biomass constraints'},length(IDs_new_variables),1);

r_sys = [r_sys;IDs_new_variables];
rNames_sys = [rNames_sys;IDs_new_variables];
ub_sys = [ub_sys; ub_new_variables];
lb_sys = [lb_sys; lb_new_variables];
obj_sys = [obj_sys; obj_new_variables];
subS_sys = [subS_sys; subS_new_variables];
grRules_sys = [grRules_sys; repmat({''}, length(IDs_new_variables),1)];

S = sparse(zeros(size(communityModel.S,1),2));

communityModel.rxns = r_sys;
communityModel.rxnNames = rNames_sys;
communityModel.ub = ub_sys;
communityModel.lb = lb_sys;
communityModel.c = obj_sys;
communityModel.subSystems = subS_sys;
communityModel.grRules_sys = grRules_sys;
communityModel.S = [communityModel.S, S];

%% add constraints
S_new_constraints = sparse(1000,size(communityModel.S,2));
n_cont_new_constraints = 0;

sense_new_contraints = repmat('E',1000,1);
new_mets = cell(1000*2,1);
lbs = cell(length(models),1);
ubs = cell(length(models),1);

for i = 1:length(models)
    lbs{i} = models{i}.lb;
    ubs{i} = models{i}.ub;
end

for i = 1:length(communityModel.rxns)
    if communityModel.ub(i)>0 && communityModel.lb(i)>=0
        communityModel.ub(i) = 99999;
        communityModel.lb(i) = 0;
    elseif communityModel.lb(i)<0 && communityModel.ub(i)<=0
        communityModel.ub(i) = 0;
        communityModel.lb(i) = -99999;
    elseif communityModel.lb(i)<0 && communityModel.ub(i)>0
        communityModel.ub(i) = 99999;
        communityModel.lb(i) = -99999;
    else
        communityModel.lb(i) = 0;
        communityModel.ub(i) = 0;
    end
end


fprintf('...adding constraints\n')
tic
for i = 1:n_microorganisms
    for k = 1:length(models{i}.rxns)
        
        if lbs{i}(k)~= 0 && lbs{i}(k)~= -99999
            %add one constraint of the form
            pos_biomass = getPosOfElementsInArray({['biomass_m' num2str(i)]},r_sys);
            pos_rxn = getPosOfElementsInArray({['model' num2str(i) '_' models{i}.rxns{k}]},r_sys);
            n_cont_new_constraints = n_cont_new_constraints+1;
            %1) Vi - lb_i*X_i >= 0
            %however it should be written as
            %1) -Vi + lb_i*X_i<= 0
            S_new_constraint = zeros(1, size(communityModel.S,2));
            S_new_constraint(pos_rxn) = 1;
            S_new_constraint(pos_biomass) = -lbs{i}(k);
            S_new_constraints(n_cont_new_constraints,:) =  S_new_constraint;
            sense_new_contraints(n_cont_new_constraints) = 'G';
            new_mets{n_cont_new_constraints} = ['new_l_cons_' num2str(i) '_' num2str(k)];
        end
        
        if ubs{i}(k)~= 0 && ubs{i}(k)~= 99999
            %add one constraint of the form
            pos_biomass = getPosOfElementsInArray({['biomass_m' num2str(i)]},r_sys);
            pos_rxn = getPosOfElementsInArray({['model' num2str(i) '_' models{i}.rxns{k}]},r_sys);
            n_cont_new_constraints = n_cont_new_constraints+1;
            %add one constraint of the form
            %2) Vi - ub_i*X_i*deltaT <= 0
            S_new_constraint = zeros(1, size(communityModel.S,2));
            S_new_constraint(pos_rxn) = 1;
            S_new_constraint(pos_biomass) = -ubs{i}(k);
            S_new_constraints(n_cont_new_constraints,:) =  S_new_constraint;
            sense_new_contraints(n_cont_new_constraints) = 'L';
            new_mets{n_cont_new_constraints} = ['new_u_cons_' num2str(i) '_' num2str(k)];
            
        end
    end
end
toc

S_new_constraints = S_new_constraints(1:n_cont_new_constraints,:);
sense_new_contraints = sense_new_contraints(1:n_cont_new_constraints,:);
new_mets = new_mets(1:n_cont_new_constraints);

communityModel.S = [communityModel.S; sparse(S_new_constraints)];
communityModel.csense = [communityModel.csense; sense_new_contraints];
communityModel.mets = [communityModel.mets; new_mets];
communityModel.metNames = [communityModel.metNames; new_mets];
communityModel.b = [communityModel.b; zeros(length(new_mets),1)];


%% add relationship between biomass and growth rate
new_mets = strcat('biomass_growth_',arrayfun(@(x) ['m_' num2str(x)], 1:n_microorganisms, 'UniformOutput', false));
S_new_constraints = zeros(n_microorganisms, size(communityModel.S,2));
sense_new_contraints = repmat('E',n_microorganisms,1);
for i = 1:n_microorganisms
    pos_biomass = getPosOfElementsInArray({['biomass_m' num2str(i)]},communityModel.rxns);
    pos_growth = getPosOfElementsInArray({['model', num2str(i),'_', biomassRxns{i}]},communityModel.rxns);
    S_new_constraints(i,pos_biomass) = -mu;
    S_new_constraints(i,pos_growth) = 1;
    sense_new_contraints(i) = 'E';
end

communityModel.S = [communityModel.S; sparse(S_new_constraints)];
communityModel.csense = [communityModel.csense; sense_new_contraints];
communityModel.mets = [communityModel.mets; new_mets'];
communityModel.metNames = [communityModel.metNames; new_mets'];
communityModel.b = [communityModel.b; zeros(length(new_mets),1)];

%% add biomass
new_mets = {'biomass_cons'};
S_new_constraints = zeros(1, size(communityModel.S,2));
sense_new_contraints = repmat('E',1,1);

for i = 1:n_microorganisms
    pos_biomass = getPosOfElementsInArray({['biomass_m' num2str(i)]},communityModel.rxns);
    S_new_constraints(pos_biomass) = 1;
end
communityModel.S = [communityModel.S; sparse(S_new_constraints)];
communityModel.csense = [communityModel.csense; sense_new_contraints];
communityModel.mets = [communityModel.mets; new_mets];
communityModel.metNames = [communityModel.metNames; new_mets];
communityModel.b = [communityModel.b; 1];

end