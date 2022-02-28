function [model, rxnsMedia, compounds, rxnsPerCompound]= setMediaFromCB(model, compounds, value)

if nargin < 3
   value = 10; 
end

all_compounds = {};
rxnsPerCompound = cell(size(compounds));
for i = 1:length(compounds)
    all_compounds = union(all_compounds, strsplit(compounds{i},';'));
    compounds_i = setdiff(strsplit(compounds{i},';'),'');
    
    if all(cellfun(@(x) ~metHasCompartment(x, {'e'}), compounds_i))
        compounds_i = addCompartmentFromMets(compounds_i, 'e', model); 
    end
    pos = cellfun(@(x) getExchangeRxnFromMetabolite(model, x), compounds_i,'UniformOutput', 0);
    if ~isempty(pos{1})
        exc = model.rxns(pos{1});
    else
        exc = {};
    end

    rxnsPerCompound{i} = exc;
end
all_compounds = setdiff(all_compounds,'');
% rxnsMedia = strcat('EX_',setdiff(rxnsMedia,''),'_e');
pos = cellfun(@(x) getExchangeRxnFromMetabolite(model, x), addCompartmentFromMets(all_compounds, 'e', model),'UniformOutput', 0);
rxnsMedia = model.rxns( cell2mat(pos(cellfun(@isempty, pos)==0))');
lowerValuesBasalMedium = -ones(size(rxnsMedia))*value;

pos_Ex = find(cellfun(@isempty, strfind(model.rxns,'EX_'))==0);
model = changeRxnBounds(model, model.rxns(pos_Ex), 0, 'l');
model = changeRxnBounds(model, model.rxns(pos_Ex(find(model.ub(pos_Ex)<0))), 0, 'u');

model = changeRxnBounds(model, rxnsMedia, lowerValuesBasalMedium, 'l');

end