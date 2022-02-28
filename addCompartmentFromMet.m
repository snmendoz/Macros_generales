function met = addCompartmentFromMet(met, compartment, modelRef)

if allMetsInCBMPYFormat(modelRef.mets)
    met = [met '_' compartment];
elseif allMetsInCOBRAFormat(modelRef.mets)
    met = [met '[' compartment ']'];
end

end