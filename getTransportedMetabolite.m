function [met,extracellular_in_left] = getTransportedMetabolite(model, rxn)

mets = findMetsFromRxns(model, rxn);
stoi = full(model.S(getPosOfElementsInArray(mets,model.mets),getPosOfElementsInArray({rxn},model.rxns)));
mets_wc = removeCompartmentFromMets(mets);

met = {};
extracellular_in_left = 0;
if length(mets)==2 &&length(unique(mets_wc))==1
    %diffusion
    met = unique(mets_wc);
elseif length(mets)==7 && all(ismember({'atp_c','adp_c','pi_c','h2o_c','h_c'},mets))
    %abc transporter
    met = unique(removeCompartmentFromMets(setdiff(mets, {'atp_c','adp_c','pi_c','h2o_c','h_c'})));
elseif length(mets)==4 && all(ismember({'pyr_c','pep_c'},mets))
    %pep transporter
    met = unique(removeCompartmentFromMets(setdiff(mets, {'pyr_c','pep_c'})));
elseif length(mets)==4 && all(ismember({'h_c','h_e'},mets))
    %symport or antiport
   met = unique(removeCompartmentFromMets(setdiff(mets, {'h_c','h_e'})));
elseif length(mets)==3 && all(ismember({'h_e'},mets))
    met = unique(removeCompartmentFromMets(setdiff(mets, {'h_e'})));
elseif length(mets)==3 && all(ismember({'h_c'},mets))
    met = unique(removeCompartmentFromMets(setdiff(mets, {'h_c'})));
end
if isempty(met); return; end;
    
if stoi(find(strcmp(mets, [met{1} '_e']))) <0
    extracellular_in_left = 1;
else
    extracellular_in_left = 0;
end

end