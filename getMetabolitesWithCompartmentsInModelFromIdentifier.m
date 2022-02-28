function mets = getMetabolitesWithCompartmentsInModelFromIdentifier(model, identifier)

mets = {};
mets_wo_comps = removeCompartmentFromMets(model.mets);
pos = find(strcmp(mets_wo_comps,identifier));
if ~isempty(pos)
    mets = model.mets(pos);
end

end