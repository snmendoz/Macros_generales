function exportModelFromCobraToJson(model,fileName)

if ~isfield(model,'grRules')
   model = creategrRulesField(model); 
end
if any(isnan(model.metCharges))
    posNaN = find(isnan(model.metCharges));
    model.metCharges(posNaN) = zeros(size(posNaN));
end
if allMetsInCOBRAFormat(model.mets)
    model = transformModelToCBMPYFormat(model);
end

fi=fopen([ fileName '.json'],'w+');

comps = getCompartmentsFromMetList(model.mets);
compartments = unique(comps);
% compartments_names =
fprintf(fi, '{\n');
%% metabolites
fprintf(fi, '"metabolites":[\n');
for i = 1:length(model.mets)-1
    fprintf(fi, '{\n');
    fprintf(fi,[ '"id":"' model.mets{i} '",\n']);
    fprintf(fi,[ '"name":"' model.metNames{i} '",\n']);
    fprintf(fi,[ '"compartment":"' comps{i} '",\n']);
    fprintf(fi,[ '"charge":' num2str(model.metCharges(i)) ',\n']);
    fprintf(fi,[ '"formula":"' model.metFormulas{i} '"\n']);
    fprintf(fi, '},\n');
end

i = length(model.mets);
fprintf(fi, '{\n');
fprintf(fi,[ '"id":"' model.mets{i} '",\n']);
fprintf(fi,[ '"name":"' model.metNames{i} '",\n']);
fprintf(fi,[ '"compartment":"' comps{i} '",\n']);
fprintf(fi,[ '"charge":' num2str(model.metCharges(i)) ',\n']);
fprintf(fi,[ '"formula":"' model.metFormulas{i} '"\n']);
fprintf(fi, '}\n');

fprintf(fi, '],\n');
%% reactions
fprintf(fi, '"reactions":[\n');
for i = 1:length(model.rxns)-1
    fprintf(fi, '{\n');
    fprintf(fi,[ '"id":"' model.rxns{i} '",\n']);
    fprintf(fi,[ '"name":"' model.rxnNames{i} '",\n']);
    pos_mets_i = find(model.S(:,i));
    mets_i = model.mets(pos_mets_i);
    coefs_i = full(model.S(pos_mets_i,i));
    fprintf(fi,[ '"metabolites":{\n']);
    for j =1:length(mets_i)-1
        if floor(coefs_i(j)) == coefs_i(j)
            fprintf(fi,[ '"' mets_i{j} '":' num2str(coefs_i(j)) '.0,\n']);
        else
            fprintf(fi,[ '"' mets_i{j} '":' num2str(coefs_i(j)) ',\n']);
        end
    end
    j = length(mets_i);
    if floor(coefs_i(j)) == coefs_i(j)
        fprintf(fi,[ '"' mets_i{j} '":' num2str(coefs_i(j)) '.0\n']);
    else
        fprintf(fi,[ '"' mets_i{j} '":' num2str(coefs_i(j)) '\n']);
    end
    fprintf(fi, '},\n');
    
    if floor(model.lb(i)) == model.lb(i)
        fprintf(fi,[ '"lower_bound":' num2str(model.lb(i)) '.0,\n']);
    else
        fprintf(fi,[ '"lower_bound":' num2str(model.lb(i)) ',\n']);
    end
    if floor(model.ub(i)) == model.ub(i)
        fprintf(fi,[ '"upper_bound":' num2str(model.ub(i)) '.0,\n']);
    else
        fprintf(fi,[ '"upper_bound":' num2str(model.ub(i)) ',\n']);
    end
    
    if i == find(model.c)
        fprintf(fi,[ '"gene_reaction_rule":"' model.grRules{i} '",\n']);
        
        objective_coeff = model.c(find(model.c));
        
        if floor(objective_coeff) == objective_coeff
            fprintf(fi,[ '"objective_coefficient":' num2str(model.c(find(model.c))) '.0\n']);
        else
            fprintf(fi,[ '"objective_coefficient":' num2str(model.c(find(model.c))) '\n']);
        end
    else
        fprintf(fi,[ '"gene_reaction_rule":"' model.grRules{i} '"\n']);
    end
    
    fprintf(fi, '},\n');
end

i = length(model.rxns);
fprintf(fi, '{\n');
fprintf(fi,[ '"id":"' model.rxns{i} '",\n']);
fprintf(fi,[ '"name":"' model.rxnNames{i} '",\n']);
pos_mets_i = find(model.S(:,i));
mets_i = model.mets(pos_mets_i);
coefs_i = full(model.S(pos_mets_i,i));
fprintf(fi,[ '"metabolites":{\n']);
for j =1:length(mets_i)-1
    if floor(coefs_i(j)) == coefs_i(j)
        fprintf(fi,[ '"' mets_i{j} '":' num2str(coefs_i(j)) '.0,\n']);
    else
        fprintf(fi,[ '"' mets_i{j} '":' num2str(coefs_i(j)) ',\n']);
    end
end
j = length(mets_i);
if floor(coefs_i(j)) == coefs_i(j)
    fprintf(fi,[ '"' mets_i{j} '":' num2str(coefs_i(j)) '.0\n']);
else
    fprintf(fi,[ '"' mets_i{j} '":' num2str(coefs_i(j)) '\n']);
end
fprintf(fi, '},\n');
if floor(model.lb(i)) == model.lb(i)
    fprintf(fi,[ '"lower_bound":' num2str(model.lb(i)) '.0,\n']);
else
    fprintf(fi,[ '"lower_bound":' num2str(model.lb(i)) ',\n']);
end
if floor(model.ub(i)) == model.ub(i)
    fprintf(fi,[ '"upper_bound":' num2str(model.ub(i)) '.0,\n']);
else
    fprintf(fi,[ '"upper_bound":' num2str(model.ub(i)) ',\n']);
end
if i == find(model.c)
    fprintf(fi,[ '"gene_reaction_rule":"' model.grRules{i} '",\n']);
    if floor(objective_coeff) == objective_coeff
        fprintf(fi,[ '"objective_coefficient":' num2str(model.c(find(model.c))) '.0\n']);
    else
        fprintf(fi,[ '"objective_coefficient":' num2str(model.c(find(model.c))) '\n']);
    end
else
    fprintf(fi,[ '"gene_reaction_rule":"' model.grRules{i} '"\n']);
end
fprintf(fi, '}\n');

fprintf(fi, '],\n');

%% genes
fprintf(fi, '"genes":[\n');
for i = 1:length(model.genes)-1
    fprintf(fi, '{\n');
    fprintf(fi,[ '"id":"' model.genes{i} '",\n']);
    fprintf(fi,[ '"name":"",\n']);
    fprintf(fi,[ '"annotation":{\n']);
    fprintf(fi,[ '"sbo":"SBO:0000243"\n']);
    fprintf(fi, '}\n');
    fprintf(fi, '},\n');
end

i = length(model.genes);
fprintf(fi, '{\n');
fprintf(fi,[ '"id":"' model.genes{i} '",\n']);
fprintf(fi,[ '"name":"",\n']);
fprintf(fi,[ '"annotation":{\n']);
fprintf(fi,[ '"sbo":"SBO:0000243"\n']);
fprintf(fi, '}\n');
fprintf(fi, '}\n');

fprintf(fi, '],\n');
%% info
fprintf(fi,[ '"id":"COBRAModel",\n']);
fprintf(fi,[ '"name":"Model Exported from COBRA Toolbox",\n']);
fprintf(fi,[ '"compartments":{\n']);
for i =1:length(compartments)-1
    fprintf(fi,[ '"' compartments{i} '":"' compartments{i} '",\n']);
end
i = length(compartments);
fprintf(fi,[ '"' compartments{i} '":"' compartments{i} '"\n']);
fprintf(fi,[ '},\n']);
fprintf(fi,[ '"version":"1"\n']);
fprintf(fi,[ '}\n']);
fclose(fi);


end