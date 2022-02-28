function [mets, pos] = getMetabolitesInRigthFromRxn(model,pos)

pos = find(model.S(:,pos)>0);
mets = model.mets(pos);

end