function [soportes, unionSoportes, posUnionSoportes] = getSoportes(model, minimos, maximos, tol)

if nargin < 4
    tol = 10^-9;
end

% para cada tramo
unionSoportes = {};

for i =1:size(minimos,1)
    minimos_i = minimos(i,:);
    maximos_i = maximos(i,:);
    soporte_i = getActivas(model, minimos_i,maximos_i, tol);
    soportes{i} = soporte_i;
    unionSoportes = union(unionSoportes, soporte_i);
end

posUnionSoportes = find(ismember(model.rxns,unionSoportes) == 1);
unionSoportes = model.rxns(posUnionSoportes);

end