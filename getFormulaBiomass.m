function getFormulaBiomass

load('iMM904.mat')
model = iMM904;
model.metCharges = model.metCharge;
pos_biomass_rxn = find(model.c);
biomass_rxn = model.rxns(pos_biomass_rxn);

pos_mets = find(model.S(:,getPosOfElementsInArray(biomass_rxn, model.rxns)));
met_info_0 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];

model = addMetabolite(model,'biomass_c');
pos_biomass_met = getPosOfElementsInArray({'biomass_c'}, model.mets);
model = addReaction(model,'growthRate','biomass_c -> ');
model.S(pos_biomass_met, pos_biomass_rxn) = 1;
getRxn_cobraFormat(model,biomass_rxn{1},1)

model = changeObjective(model, 'growthRate');
fba = optimizeCbModel(model);

[model2, elements, rCoefs] = inferMetFormula(model,biomass_rxn{1},1);
pos_mets = find(model2.S(:,getPosOfElementsInArray(biomass_rxn, model2.rxns)));
met_info_2 = [model2.mets(pos_mets) model2.metFormulas(pos_mets)];

coef_biomasss = [elements' num2cell(rCoefs*1000) num2cell(rCoefs) num2cell(rCoefs/rCoefs(1))];
formula1 = ['C' num2str(coef_biomasss{1,3}) 'H' num2str(coef_biomasss{2,3}) 'O' num2str(coef_biomasss{3,3}) 'N' num2str(coef_biomasss{4,3}) 'P' num2str(coef_biomasss{5,3}) 'S' num2str(coef_biomasss{6,3})];
formula2 = ['C' num2str(coef_biomasss{1,2}) 'H' num2str(coef_biomasss{2,2}) 'O' num2str(coef_biomasss{3,2}) 'N' num2str(coef_biomasss{4,2}) 'P' num2str(coef_biomasss{5,2}) 'S' num2str(coef_biomasss{6,2})];
MW = calculateMolecularWeight_fromFormula(['C' num2str(coef_biomasss{1,2}) 'H' num2str(coef_biomasss{2,2}) 'O' num2str(coef_biomasss{3,2}) 'N' num2str(coef_biomasss{4,2}) 'P' num2str(coef_biomasss{5,2}) 'S' num2str(coef_biomasss{6,2})])
% MW = 963520 g/mol
%% biomass 'C36H66O24N6': 
%C: 1
%H 1.8333
%O: 0.6667
%N: 0.1667

% C36H64O18N5, con modelo yeast8
end