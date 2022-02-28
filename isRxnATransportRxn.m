function tf = isRxnATransportRxn(model, pos)

tf = arrayfun(@(x) length(unique(getCompartmentsFromMetList(model.mets(find(model.S(:,pos(x)))))))>1, 1:length(pos));

end