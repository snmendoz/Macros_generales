function activas = getActivas(model, minimos, maximos, tol)

if nargin < 4
    tol = 10^-9;
end

posActivas = union(find(minimos > tol), find(maximos < -tol));
activas = model.rxns(posActivas);

end