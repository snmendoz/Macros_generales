function essentials = singleOmission(model, media, rxns,threshold_percentage)

if nargin<2;  media = {}; end
if nargin<3 || isempty(rxns);  rxns = model.rxns(findExcRxns(model)); end
if nargin<4 || isempty(threshold_percentage);  threshold_percentage = 0.1; end;

fba = optimizeCbModel(model);
threshold_absolute = fba.f*threshold_percentage;

if ~isempty(media)
    
end

essentials = cell(size(rxns));
n_essentials = 0;

for i = 1:length(rxns)
    modelAux = changeRxnBounds(model, rxns(i),0,'l');
    if modelAux.ub(find(strcmp(modelAux.rxns,rxns(i))))<modelAux.lb(find(strcmp(modelAux.rxns,rxns(i))))
        modelAux.ub(find(strcmp(modelAux.rxns,rxns(i)))) = 10; 
    end
    fbaAux = optimizeCbModel(modelAux);
    if isempty(fbaAux.x) || fbaAux.f<=threshold_absolute
        n_essentials = n_essentials+1;
        essentials{n_essentials} = rxns{i};
    end
end

essentials = essentials(1:n_essentials);

end