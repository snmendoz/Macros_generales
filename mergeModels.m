function [mergedModel, exchangeRxns] = mergeModels(models, sharedMetaboliteIDs, varargin)

parser = inputParser();
parser.addRequired('models', @iscell)
parser.addRequired('sharedMetaboliteIDs', @iscell)
parser.addParameter('compartmentSep', {'[',']'}, @iscell)
defaultPrefixes = {};
for i = 1:length(models)
    defaultPrefixes{i} = ['model' num2str(i)];
end
parser.addParameter('prefixes', defaultPrefixes, @(x) iscell(x) && length(models) == length(x))
parser.addParameter('description', 'community model' , @ischar)
parser.addParameter('indMets', {{},{}}, @iscell)
parser.addParameter('indLb', {[],[]}, @iscell)
parser.addParameter('indUb', {[],[]}, @iscell)
parser.addParameter('lbPool', -1000*ones(size(sharedMetaboliteIDs)), @isnumeric)
parser.addParameter('ubPool', 1000*ones(size(sharedMetaboliteIDs)), @isnumeric)

parser.parse(models, sharedMetaboliteIDs, varargin{:})
models = parser.Results.models;
sharedMetaboliteIDs = parser.Results.sharedMetaboliteIDs;
compartmentSep = parser.Results.compartmentSep;
prefixes = parser.Results.prefixes;
description = parser.Results.description;
indMets = parser.Results.indMets;
indLb = parser.Results.indLb;
indUb = parser.Results.indUb;
lbPool = parser.Results.lbPool;
ubPool = parser.Results.ubPool;

for i = 1:length(models)
    if isempty(indUb{i}) && ~isempty(indMets{i})
        indUb{i} = 1000*ones(size(indMets{i}));
    end
    if isempty(indLb{i}) && ~isempty(indMets{i})
        indLb{i} = -1000*ones(size(indMets{i}));
    end
end

n_mets = [];
n_rxns = [];
for i = 1:length(models)
    [n, m] = size(models{i}.S);
    n_mets = [n_mets; n];
    n_rxns = [n_rxns; m];
end

n_shared = length(sharedMetaboliteIDs);
for i = 1:length(models)
    n_onlyInModel(i) = length(indMets{i});
end
% create metabolites in the medium
mets_pool = strcat(sharedMetaboliteIDs, compartmentSep(1), 'pool', compartmentSep(2));
for i = 1:length(models)
    mets_pool = [mets_pool; strcat(indMets{i}, compartmentSep(1), 'pool', compartmentSep(2))];
end
% create exchange reactions for each metabolite
rxns_pool = strcat('EX_',  sharedMetaboliteIDs, '_pool');
for i = 1:length(models)
    rxns_pool = [rxns_pool; strcat('EX_', indMets{i}, '_pool')];
end
% define lower and upper bounds for exchange rxns
lb_pool = lbPool;
ub_pool = ubPool;

for i = 1:length(models)
    lb_pool = [lb_pool; indLb{i}];
    ub_pool = [ub_pool; indUb{i}];
end

% matrices for new metabolites and reactions
pool = zeros(n_shared+sum(n_onlyInModel),sum(n_rxns));
pool_ex = -1*eye(n_shared+sum(n_onlyInModel));
exchangeRxns = {};
for i = 1:length(models)
    posMets_i = cell2mat(arrayfun(@(x)find(strcmp(x, models{i}.mets)), strcat(sharedMetaboliteIDs, compartmentSep(1), 'e', compartmentSep(2) ), 'UniformOutput', false))';
    for j = 1:length(posMets_i)
        pos_ex = intersect(find(sum(models{i}.S~=0,1)==1),find(models{i}.S(posMets_i(j),:)));        
        if isempty(pos_ex); continue; end;
        if i > 1
            newPos_ex = pos_ex + sum(n_rxns(1:i-1));
        else
            newPos_ex = pos_ex;
        end
        exchangeRxns = [exchangeRxns; {[prefixes{i} '_' models{i}.rxns{pos_ex}]}];
        if models{i}.S(posMets_i(j), pos_ex) < 0
            pool(j, newPos_ex) = 1;
        else
            pool(j, newPos_ex) = -1;
        end
    end
end

posMet = n_shared;
for i = 1:length(indMets)
    mets_i = indMets{i};
    for j = 1:length(mets_i)
        posMet = posMet+1;
                
        pos_met = find(strcmp(models{i}.mets, strcat(mets_i{j}, compartmentSep(1), 'e', compartmentSep(2) )));
        
        pos_ex = intersect(find(sum(models{i}.S~=0,1)==1),find(models{i}.S(pos_met,:)));
        
        if i > 1
            newPos_ex = pos_ex + sum(n_rxns(1:i-1));
        else
            newPos_ex = pos_ex;
        end
        exchangeRxns = [exchangeRxns; {[prefixes{i} '_' models{i}.rxns{pos_ex}]}];
        if models{i}.S(pos_met, pos_ex) < 0
            pool(posMet, newPos_ex) = 1;
        else
            pool(posMet, newPos_ex) = -1;
        end        
    end
end

% System Stoichiometric Matrix
S_sys = full(models{1}.S);
for i = 2:length(models)
    S_sys = blkdiag(S_sys, full(models{i}.S));
end
S_sys = [S_sys, zeros(sum(n_mets),n_shared+sum(n_onlyInModel)); ...
         pool , pool_ex];

% system metabolites
m_sys = strcat(prefixes{1},'_',models{1}.mets);
mNames_sys = strcat(prefixes{1},'_',models{1}.metNames);
for i = 2:length(models)
    m_sys = [m_sys; strcat(prefixes{i},'_',models{i}.mets)];
    mNames_sys = [mNames_sys; strcat(prefixes{i},'_',models{i}.metNames)];
end
m_sys = [m_sys; mets_pool];
mNames_sys = [mNames_sys; mets_pool];

% system reactions
r_sys = strcat(prefixes{1},'_',models{1}.rxns);
rNames_sys = strcat(prefixes{1},'_',models{1}.rxnNames);
ub_sys = models{1}.ub;
lb_sys = models{1}.lb;
obj_sys = models{1}.c;
subS_sys = models{1}.subSystems;
for i = 2:length(models)
    r_sys = [r_sys; strcat(prefixes{i},'_',models{i}.rxns)];
    rNames_sys = [rNames_sys; strcat(prefixes{i},'_',models{i}.rxnNames)];
    ub_sys = [ub_sys; models{i}.ub];
    lb_sys = [lb_sys; models{i}.lb];
    obj_sys = [obj_sys; models{i}.c];
    subS_sys = [subS_sys; models{i}.subSystems];
end
r_sys = [r_sys; rxns_pool];
rNames_sys = [rNames_sys; rxns_pool];
ub_sys = [ub_sys; ub_pool];
lb_sys = [lb_sys; lb_pool];
obj_sys = [obj_sys; zeros(n_shared+sum(n_onlyInModel),1)];
subS_sys = [subS_sys; cell(n_shared+sum(n_onlyInModel),1)];

pos = find(cellfun(@isempty, subS_sys));
for i=1:length(pos)
    subS_sys{pos(i)} = '';
end

subS_sys = regexprep(subS_sys, {'transport', 'Nucleotide Metabolism'}, {'Transport', 'Nucleotides Metabolism'});

% create model
mergedModel = struct('S', S_sys, 'mets', {m_sys}, 'metNames', {mNames_sys},...
    'subSystems', {subS_sys}, ...
    'rxns', {r_sys}, 'rxnNames', {rNames_sys}, 'b', zeros(size(m_sys)), ...
    'ub', ub_sys, 'lb', lb_sys, 'c', obj_sys, 'description', description);

end