function [mets, pos] = getMetabolitesInRightFromRxn(model,pos)

pos = find(model.S(:,pos)>0);
mets = model.mets;

end