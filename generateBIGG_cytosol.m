function names = generateBIGG_cytosol(database)
load('models_refined.mat')
load(database);
removed=cell(10000,1);
removed_p=cell(10000,1);
cont=1;
cont_p=1;
biggGeneral=bigg;

for i=1:length(biggGeneral.rxns)
    mets=biggGeneral.mets(find(biggGeneral.S(:,i)));
    if length(find(cellfun(@isempty,strfind(mets,'[c]'))==0)) ~=length(mets)
        removed(cont)=biggGeneral.rxns(i);
        cont=cont+1;
    end
    if length(find(cellfun(@isempty,strfind(mets,'[c]'))==0)) + length(find(cellfun(@isempty,strfind(mets,'[p]'))==0)) ~=length(mets)
        removed_p(cont_p)=biggGeneral.rxns(i);
        cont_p=cont_p+1;
    end
end

model=biggGeneral;
removed=removed(1:cont-1);
biggGeneral=removeRxns(model,removed);
bigg = biggGeneral;
save([database '_cytosol'],'bigg');
% cd C:\Users\notebook\Desktop\Trabajo\cobra_antiguo\io
% writeCbModel2(bigg,'sbml','bigg_80_fixed_forMeneco_cytosol');
% cd('2C:\Users\notebook\Desktop\Potencial Redox\Oenococcus oeni\Macros\BIGG');

model=biggGeneral;
removed_p=removed_p(1:cont_p-1);
biggGeneral=removeRxns(model,removed_p);
bigg = biggGeneral;
save([database '_cytosol_periplasm'],'bigg');

names{1} = [database '_cytosol'];
names{2} = [database '_cytosol_periplasm'];

% cd C:\Users\notebook\Desktop\Trabajo\cobra_antiguo\io
% writeCbModel2(bigg,'sbml','bigg_80_fixed_forMeneco_cytosol_periplasm');
% cd('2C:\Users\notebook\Desktop\Potencial Redox\Oenococcus oeni\Macros\BIGG');

%create database for bacteria, in general
%identify models_refined which are bacteria
bact_bool = ones(size(models_refined));
bact_grampos_bool = zeros(size(models_refined));
bact_gramneg_bool = zeros(size(models_refined));

rxns_bacteria = {};
rxns_grampos = {};
rxns_gramneg = {};

for i = 1:length(models_refined)
    model_i = models_refined{i};
    [tok,rem] = strtok(model_i.mets,'[');
    compartments = unique(regexprep(rem,{'[',']'},{'',''}));
    if ~isempty(setdiff(compartments,{'c','p','e'}))
        bact_bool(i)=0;
    end
    if length(intersect({'c','e'}, compartments)) == length(compartments) && length(compartments) == 2
        bact_grampos_bool(i)=1;
    end
    if length(intersect({'c','e','p'}, compartments)) == length(compartments) && length(compartments) == 3
        bact_gramneg_bool(i)=1;
    end
    if bact_bool(i);
        rxns_bacteria = union(rxns_bacteria,model_i.rxns);
    end
    if bact_grampos_bool(i);
        rxns_grampos = union(rxns_grampos,model_i.rxns);
    end
    if bact_gramneg_bool(i);
        rxns_gramneg = union(rxns_gramneg,model_i.rxns);
    end
end
biggBacteria = removeRxns(biggGeneral, setdiff(biggGeneral.rxns,rxns_bacteria));
bigg = biggBacteria;
save([database '_bacteria'],'bigg');
biggGramPos = removeRxns(biggGeneral, setdiff(biggGeneral.rxns,rxns_grampos));
bigg = biggGramPos;
save([database '_grampos'],'bigg');
biggGramNeg = removeRxns(biggGeneral, setdiff(biggGeneral.rxns,rxns_gramneg));
bigg = biggGramNeg;
save([database '_gramneg'],'bigg');

end