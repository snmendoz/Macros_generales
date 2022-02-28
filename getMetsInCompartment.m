function mets = getMetsInCompartment(model, compartment)

comps = getCompartmentsFromMetList(model.mets);
mets = model.mets(strcmp(comps,compartment));

end