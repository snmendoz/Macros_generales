function model = addRateRatioToModel(model, rxn1, rxn2, ratioName, ratio, sameSense)

%rxn1/rxn2 = ratio

if sameSense
   ratio = -ratio; 
end

model = addMetabolite(model, ratioName);
model.S(find(strcmp(model.mets, ratioName)), find(strcmp(model.rxns,rxn1)) )=1;
model.S(find(strcmp(model.mets, ratioName)), find(strcmp(model.rxns,rxn2)) )=ratio;

end