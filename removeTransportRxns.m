function model = removeTransportRxns(model, comp1, comp2)

model2 = model;
if isempty(strfind(model2.mets{1}, '['))
    model2.mets = regexprep(model2.mets, '_(.)$', '[$1]');  
end
if isempty(strfind(comp1, '['))
    comp1 = regexprep(comp1,'_(.)$','[$1]');
end
if isempty(strfind(comp2, '['))
    comp2 = regexprep(comp2,'_(.)$','[$1]');
end
    
[compartmentReactions1] = findRxnFromCompartment(model2, comp1);
rxns1 = compartmentReactions1(:,1);
[compartmentReactions2] = findRxnFromCompartment(model2, comp2);
rxns2 = compartmentReactions2(:,1);

model = removeRxns(model, intersect(rxns1,rxns2));

end