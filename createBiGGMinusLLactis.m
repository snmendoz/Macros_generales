function createBiGGMinusLLactis(bigg, iNF517)

if nargin<1
    load('D:\Dropbox\Databases\BIGG\bigg_85.mat')
end
if nargin<2
    load('D:\Dropbox\Databases\BIGG\iNF517.mat');
end
% iNF517 is the 61
load('D:\Dropbox\Databases\BIGG\Presence_rxns_bigg_85.mat')
% which reactions are only in iNF517
pos_rxnsIniNF517 = find(sum(Presence_rxns(:,61)==1,2)==1);
pos_rxnsOnlyInOneModel = find(sum(Presence_rxns==1,2)==1);
pos_rxnsOnlyInINF517 = intersect(pos_rxnsIniNF517,pos_rxnsOnlyInOneModel);

bigg_minus_iNF517 = removeRxns(bigg, bigg.rxns(pos_rxnsOnlyInINF517));
save('D:\Dropbox\Databases\BIGG\bigg_85_minus_iNF517.mat', 'bigg_minus_iNF517');

end