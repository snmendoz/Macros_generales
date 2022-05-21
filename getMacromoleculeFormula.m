function getMacromoleculeFormula

cd('C:\Users\snmen\Dropbox\Research_Projects\YeastGEM')
load('Yeast833_bigg_genesFixed')
%% biomass

getRxn_cobraFormat(model,'r_4041',1)
getRxn_cobraFormat(model,'r_4041',0)

pos_mets = find(model.S(:,getPosOfElementsInArray({'r_4041'}, model.rxns)));
met_info_0 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];

%% protein reaction: r_4047 (protein: s_3717[c])
model.rxns(find(model.S(getPosOfElementsInArray({'s_3717[c]'},model.mets),:)))
getRxn_cobraFormat(model,'r_4047',1)
% '0.52701 Ala-tRNA(Ala) [cytoplasm] + 0.18459 Arg-tRNA(Arg) [cytoplasm] + 0.11682 Asn-tRNA(Asn) [cytoplasm] + 0.34173 Asp-tRNA(Asp) [cytoplasm] + 0.0075813 Cys-tRNA(Cys) [cytoplasm] + 0.12107 Gln-tRNA(Gln) [cytoplasm] + 0.34667 Glu-tRNA(Glu) [cytoplasm] + 0.33358 Gly-tRNA(Gly) [cytoplasm] + 0.076157 His-tRNA(His) [cytoplasm] + 0.22135 Ile-tRNA(Ile) [cytoplasm] + 0.34047 Leu-tRNA(Leu) [cytoplasm] + 0.32875 Lys-tRNA(Lys) [cytoplasm] + 0.058238 Met-tRNA(Met) [cytoplasm] + 0.15381 Phe-tRNA(Phe) [cytoplasm] + 0.18919 Pro-tRNA(Pro) [cytoplasm] + 0.21296 Ser-tRNA(Ser) [cytoplasm] + 0.21986 Thr-tRNA(Thr) [cytoplasm] + 0.032622 Trp-tRNA(Trp) [cytoplasm] + 0.11716 Tyr-tRNA(Tyr) [cytoplasm] + 0.30394 Val-tRNA(Val) [cytoplasm] -> 0.52701 tRNA(Ala) [cytoplasm] + 0.18459 tRNA(Arg) [cytoplasm] + 0.11682 tRNA(Asn) [cytoplasm] + 0.34173 tRNA(Asp) [cytoplasm] + 0.0075813 tRNA(Cys) [cytoplasm] + 0.12107 tRNA(Gln) [cytoplasm] + 0.34667 tRNA(Glu) [cytoplasm] + 0.33358 tRNA(Gly) [cytoplasm] + 0.076157 tRNA(His) [cytoplasm] + 0.22135 tRNA(Ile) [cytoplasm] + 0.34047 tRNA(Leu) [cytoplasm] + 0.32875 tRNA(Lys) [cytoplasm] + 0.058238 tRNA(Met) [cytoplasm] + 0.15381 tRNA(Phe) [cytoplasm] + 0.18919 tRNA(Pro) [cytoplasm] + 0.21296 tRNA(Ser) [cytoplasm] + 0.21986 tRNA(Thr) [cytoplasm] + 0.032622 tRNA(Trp) [cytoplasm] + 0.11716 tRNA(Tyr) [cytoplasm] + 0.30394 tRNA(Val) [cytoplasm] + 1 protein [cytoplasm]'
pos_mets = find(model.S(:,getPosOfElementsInArray({'r_4047'}, model.rxns)));
met_info_1 = [model.mets(pos_mets) model.metFormulas(pos_mets)];

model = inferMetFormula(model,'r_4047',0);
pos_mets = find(model.S(:,getPosOfElementsInArray({'r_4047'}, model.rxns)));
met_info_2 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];

%'C20H36O6N6' es igual a CH1.8O0.3N0.3

%% carbohidrate (s_3718[c]). Carbohydrate rxn
model.rxns(find(model.S(getPosOfElementsInArray({'s_3718[c]'},model.mets),:)))
getRxn_cobraFormat(model,'r_4048',1)
pos_mets = find(model.S(:,getPosOfElementsInArray({'r_4048'}, model.rxns)));
met_info_3 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];

model.metFormulas{getPosOfElementsInArray({'mannan[c]'}, model.mets)} = 'C6H10O5';
model.metFormulas{getPosOfElementsInArray({'glycogen[c]'}, model.mets)} = 'C6H10O5';
met_info_3 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];
model = inferMetFormula(model,'r_4048',0);
pos_mets = find(model.S(:,getPosOfElementsInArray({'r_4048'}, model.rxns)));
met_info_4 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];

%% RNA(s_3719[c])
model.rxns(find(model.S(getPosOfElementsInArray({'s_3719[c]'},model.mets),:)))
rxn = 'r_4049';
getRxn_cobraFormat(model,rxn,1)
pos_mets = find(model.S(:,getPosOfElementsInArray({rxn}, model.rxns)));
met_info_3 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];

model = inferMetFormula(model,rxn,0);
pos_mets = find(model.S(:,getPosOfElementsInArray({rxn}, model.rxns)));
met_info_4 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];

%% DNA reaction: (s_3720[c])

model.rxns(find(model.S(getPosOfElementsInArray({'s_3720[c]'},model.mets),:)))
rxn = 'r_4050';
getRxn_cobraFormat(model,rxn,1)
pos_mets = find(model.S(:,getPosOfElementsInArray({rxn}, model.rxns)));
met_info_3 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];

model = inferMetFormula(model,rxn,0);
pos_mets = find(model.S(:,getPosOfElementsInArray({rxn}, model.rxns)));
met_info_4 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];


%% cofactor: (s_3720[c])

model.rxns(find(model.S(getPosOfElementsInArray({'s_4205[c]'},model.mets),:)))
rxn = 'r_4598';
getRxn_cobraFormat(model,rxn,1)
pos_mets = find(model.S(:,getPosOfElementsInArray({rxn}, model.rxns)));
met_info_3 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];

model = inferMetFormula(model,rxn,0);
pos_mets = find(model.S(:,getPosOfElementsInArray({rxn}, model.rxns)));
met_info_4 = [model.mets(pos_mets) model.metFormulas(pos_mets)];

%% ion: (s_42060[c])
model.S(getPosOfElementsInArray({'s_4206[c]'},model.mets),getPosOfElementsInArray({'r_4041'}, model.rxns)) = 0;

% model.rxns(find(model.S(getPosOfElementsInArray({'s_4206[c]'},model.mets),:)))
% rxn = 'r_4599';
% getRxn_cobraFormat(model,rxn,1)
% pos_mets = find(model.S(:,getPosOfElementsInArray({rxn}, model.rxns)));
% met_info_3 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];
% 
% model = inferMetFormula(model,rxn,0);
% pos_mets = find(model.S(:,getPosOfElementsInArray({rxn}, model.rxns)));
% met_info_4 = [model.mets(pos_mets) model.metFormulas(pos_mets)];

%% biomass

% getRxn_cobraFormat(model,'r_4041',1)
% getRxn_cobraFormat(model,'r_4041',0)
% 
% pos_mets = find(model.S(:,getPosOfElementsInArray({'r_4041'}, model.rxns)));
% met_info_0 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];
% 
% model.metFormulas{getPosOfElementsInArray({'s_0450[c]'},model.mets)} = 'C35.942H65.5731O23.6762N5.5965P0.1977S0.0773';
% 
% rxn = 'r_4041';
% model = inferMetFormula(model,rxn,0);
% pos_mets = find(model.S(:,getPosOfElementsInArray({rxn}, model.rxns)));
% met_info_4 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];

%% lipid (s_1096[c])
model.rxns(find(model.S(getPosOfElementsInArray({'s_1096[c]'},model.mets),:)))
getRxn_cobraFormat(model,'r_2108',1)
getRxn_cobraFormat(model,'r_2108',0)

%% lipid chain
model.rxns(find(model.S(getPosOfElementsInArray({'s_3747[c]'},model.mets),:)))
getRxn_cobraFormat(model,'r_4064',1)
getRxn_cobraFormat(model,'r_4065',1)
rxn = 'r_4065';
pos_mets = find(model.S(:,getPosOfElementsInArray({'r_4065'}, model.rxns)));
met_info_0 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];
model = inferMetFormula(model,rxn,0);
pos_mets = find(model.S(:,getPosOfElementsInArray({'r_4065'}, model.rxns)));
met_info_2 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];

%% fatty acid backbone
model.rxns(find(model.S(getPosOfElementsInArray({'s_0694[c]'},model.mets),:)))
getRxn_cobraFormat(model,model.rxns(find(model.S(getPosOfElementsInArray({'s_0694[c]'},model.mets),:))),1)

rxn = 'r_3977';
pos_mets = find(model.S(:,getPosOfElementsInArray({rxn}, model.rxns)));
met_info_0 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];

model = inferMetFormula(model,rxn,0);
pos_mets = find(model.S(:,getPosOfElementsInArray({rxn}, model.rxns)));
met_info_2 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];

%% lipid backbone
model.rxns(find(model.S(getPosOfElementsInArray({'s_3746[c]'},model.mets),:)))
getRxn_cobraFormat(model,'r_4063',1)
rxn = 'r_4063';
pos_mets = find(model.S(:,getPosOfElementsInArray({rxn}, model.rxns)));
met_info_0 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];
model = inferMetFormula(model,rxn,0);
pos_mets = find(model.S(:,getPosOfElementsInArray({rxn}, model.rxns)));
met_info_2 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];

%%
model.rxns(find(model.S(getPosOfElementsInArray({'s_1096[c]'},model.mets),:)))
rxn = 'r_2108';
pos_mets = find(model.S(:,getPosOfElementsInArray({rxn}, model.rxns)));
met_info_3 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];
model = inferMetFormula(model,rxn,0);
pos_mets = find(model.S(:,getPosOfElementsInArray({rxn}, model.rxns)));
met_info_2 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];

% checkMassChargeBalance(model,-1)
%% biomass

model.rxns(find(model.S(getPosOfElementsInArray({'s_3720[c]'},model.mets),:)))
rxn ='r_4041';
getRxn_cobraFormat(model,rxn,1)
pos_mets = find(model.S(:,getPosOfElementsInArray({rxn}, model.rxns)));
met_info_3 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];

model = inferMetFormula(model,rxn,0);
pos_mets = find(model.S(:,getPosOfElementsInArray({rxn}, model.rxns)));
met_info_4 = [model.mets(pos_mets) model.metNames(pos_mets) model.metFormulas(pos_mets)];

model.metFormulas{getPosOfElementsInArray({'s_0450[c]'},model.mets)}

% 'r_3975': C36H64O18N5
% 'r_3976': C36H64O18N5
% 'r_3977': C36H64O18N5

end