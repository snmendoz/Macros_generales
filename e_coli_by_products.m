% fame
load('D:\Dropbox\Databases\BIGG\iJO1366.mat')
model = iJO1366;
model.csense = model.csense'
model = changeRxnBounds(model, 'EX_glc__D_e',-10,'l');
printConstraints(model, -1000,1000)
fba = optimizeCbModel(model);

model2 = changeObjective(model, 'EX_etoh_e');
model2 = changeRxnBounds(model2, model.rxns(model.c==1),fba.f,'l');
max_ethanol = optimizeCbModel(model2);
model2 = changeObjective(model, 'EX_ac_e');
model2 = changeRxnBounds(model2, model.rxns(model.c==1),fba.f,'l');
max_acetate= optimizeCbModel(model2);
model2 = changeObjective(model, 'EX_for_e');
model2 = changeRxnBounds(model2, model.rxns(model.c==1),fba.f,'l');
max_formate= optimizeCbModel(model2);
model2 = changeObjective(model, 'EX_lac__D_e');
model2 = changeRxnBounds(model2, model.rxns(model.c==1),fba.f,'l');
max_lactate= optimizeCbModel(model2);

disp('')