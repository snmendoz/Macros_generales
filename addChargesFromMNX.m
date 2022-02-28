function model = addChargesFromMNX(model, mets, language, metOtherIDs, metMNXCharges, metMNXFormulas)

metswc = removeCompartmentFromMets(mets);
for i = 1:length(mets)
    posInMNX = getPosOfElementsInArray({[language ':' metswc{i}]}, metOtherIDs);
    posInModel = getPosOfElementsInArray(mets(i), model.mets);
    if ~isempty(posInMNX) && ~isempty(model.metFormulas{posInModel}) && ~isempty(metMNXFormulas{posInMNX})
        if sameFormula(model.metFormulas{posInModel}, metMNXFormulas{posInMNX}) 
            model.metCharges(posInModel) = metMNXCharges(posInMNX);
        else
            if isFormulaDefined(metMNXFormulas{posInMNX})
                charge = calculateChargeFromReference(model.metFormulas{posInModel}, metMNXFormulas{posInMNX}, metMNXCharges(posInMNX));
                model.metCharges(posInModel) = charge;
            end            
        end
    end
end

end