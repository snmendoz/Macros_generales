function has = rxnHasCompartment(model, rxn, listOfCompartments)

mets = getMetaboliteIDsFromRxns(model, {rxn});
comps = getCompartmentsFromMetList(mets{1});
has = 0;
if length(intersect(comps, listOfCompartments)) == length(listOfCompartments)
    has = 1;
end
end