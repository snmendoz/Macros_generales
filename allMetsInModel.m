function true = allMetsInModel(model, mets, relaxCompartment)

if nargin <3
    relaxCompartment = 0;
end

true = 0;
if all(cellfun(@(x) metInModel(model, x, relaxCompartment), mets))
    true = 1;
end

end