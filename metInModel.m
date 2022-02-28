function true = metInModel(model, met, relaxCompartment)

if nargin <3
    relaxCompartment = 0;
end
true = 0;
if ~relaxCompartment && ismember(met, model.mets)
    true = 1;
elseif relaxCompartment && ismember(removeCompartmentFromMets({met}), removeCompartmentFromMets(model.mets))
    true = 1;
end

end