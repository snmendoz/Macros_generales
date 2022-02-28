function [scores, sets] = run_test_biolog(model, sets, insilico_cutoff)

TP_all = 0;
TN_all = 0;
FP_all = 0;
FN_all = 0;
for i = 1:length(sets)
    
    exp_results = sets{i}.exp_results;
    media = sets{i}.media;
    media_lb = sets{i}.media_lb;
    sources = sets{i}.sources;
    sources_lb = sets{i}.sources_lb;
    insilico_results = zeros(size(exp_results));
    insilico_growth = zeros(size(exp_results));
    
    
    for j = 1:length(exp_results)
        ex = model.rxns(intersect(find(model.lb<0),strmatch('EX_',model.rxns)));
        model = changeRxnBounds(model, ex, 0, 'l');
        model = changeRxnBounds(model, media, media_lb, 'l');
        model = changeRxnBounds(model, sources, 0, 'l');
        model = changeRxnBounds(model, sources{j}, sources_lb(j), 'l');
        fba_j = optimizeCbModel(model);
        insilico_growth(j) = fba_j.f;   
        if fba_j.f > insilico_cutoff
            insilico_results(j) = 1;
        end
    end
    
    pos_TP =  intersect(find(exp_results), find(insilico_results));
    pos_TN = intersect(find(exp_results==0), find(insilico_results==0));
    pos_FP = intersect(find(exp_results==0), find(insilico_results));
    pos_FN = intersect(find(exp_results), find(insilico_results==0));
        
    TP = length(pos_TP);
    TN = length(pos_TN);
    FP = length(pos_FP);
    FN = length(pos_FN);
    
    TP_all = TP_all + TP;
    TN_all = TN_all + TN;
    FP_all = FP_all + FP;
    FN_all = FN_all + FN;
    
    % Paso 7: determinar parametros de calidad del medelo e imprimir resultados
    SEN = TP / (TP + FN);
    ESP = TN / (TN + FP);
    PRE = TP / (TP + FP);
    NPV = TN / (TN + FN);
    ACC = (TP + TN) / (TN + TP + FN + FP);
    FSCORE = (2 * (PRE * SEN))/ (PRE + SEN);
    
    prediction = cell(size(sources));
    prediction(pos_TP) = repmat({'TP'},length(pos_TP),1);
    prediction(pos_TN) = repmat({'TN'},length(pos_TN),1);
    prediction(pos_FP) = repmat({'FP'},length(pos_FP),1);
    prediction(pos_FN) = repmat({'FN'},length(pos_FN),1);
    
    sets{i}.prediction = prediction;
    sets{i}.insilico_results = insilico_results;
    sets{i}.insilico_growth = insilico_growth;
    sets{i}.posTP = pos_TP;
    sets{i}.posTN = pos_TN;
    sets{i}.posFP = pos_FP;
    sets{i}.posFN = pos_FN;
    sets{i}.TP = TP;
    sets{i}.TN = TN;
    sets{i}.FP = FP;
    sets{i}.FN = FN;
    sets{i}.SEN = SEN;
    sets{i}.ESP = ESP;
    sets{i}.PRE = PRE;
    sets{i}.NPV = NPV;
    sets{i}.ACC = ACC;
    sets{i}.FSCORE = FSCORE;
    
end

SEN = TP_all / (TP_all + FN_all);
ESP = TN_all / (TN_all + FP_all);
PRE = TP_all / (TP_all + FP_all);
NPV = TN_all / (TN_all + FN_all);
ACC = (TP_all + TN_all) / (TN_all + TP_all + FN_all + FP_all);
FSCORE = (2 * (PRE * SEN))/ (PRE + SEN);

scores.TP = TP_all;
scores.TN = TN_all;
scores.FP = FP_all;
scores.FN = FN_all;
scores.SEN = SEN;
scores.ESP = ESP;
scores.PRE = PRE;
scores.NPV = NPV;
scores.ACC = ACC;
scores.FSCORE = FSCORE;

end