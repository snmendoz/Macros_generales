function [activeConstraints, notActiveConstraints, posActiveConstraints, posNotActiveConstraints] = findActiveConstraints(model,solution, tol)
if nargin <3
    tol = 1e-9;
end
posActiveConstraints = find(abs(model.S*solution - model.b)<tol);
posNotActiveConstraints = setdiff(1:length(model.mets),posActiveConstraints);
activeConstraints = model.mets(posActiveConstraints);
notActiveConstraints = model.mets(posNotActiveConstraints);

end