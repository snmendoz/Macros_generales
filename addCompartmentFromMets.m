function mets = addCompartmentFromMets(mets, compartment, modelRef)

if allMetsInCBMPYFormat(modelRef.mets)
    mets = strcat(mets, '_', compartment);
elseif allMetsInCOBRAFormat(modelRef.mets)
    mets = strcat(mets, '[', compartment,']');
end

end




