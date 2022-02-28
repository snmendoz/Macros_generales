function [n, nt, p, metsInModel, metsNotInModel] = calculateMetsInModelPerReaction(model, modelRef)

% %slowest
% tic
% n = arrayfun(@(x) length(find(ismember(model.mets(model.S(:,x)~=0), modelRef.mets)==1)), 1:length(model.rxns));
% toc

%fastest
n = zeros(size(model.rxns));
nt = zeros(size(model.rxns));
p = zeros(size(model.rxns));
metsInModel = cell(size(model.rxns));
metsNotInModel = cell(size(model.rxns));

for i = 1:length(model.rxns)    
    mets = model.mets(model.S(:,i)~=0);
    metsInModel_i = mets(ismember(mets, modelRef.mets)==1);
    metsNotInModel_i = mets(ismember(mets, modelRef.mets)==0);
    n(i) = length(metsInModel_i);
    nt(i) = length(mets);
    p(i) = 100*n(i)/nt(i);
    metsInModel{i} = metsInModel_i;
    metsNotInModel{i} = metsNotInModel_i;
end

% %third fastest
% tic
% n = arrayfun(@(x) length(intersect(model.mets(model.S(:,x)~=0), modelRef.mets)), 1:length(model.rxns));
% toc
% 
% %second fastest
% tic
% n = zeros(size(model.rxns));
% for i = 1:length(model.rxns)    
%     n(i) = length(intersect(model.mets(model.S(:,1)~=0), modelRef.mets));
% end
% tocs


end