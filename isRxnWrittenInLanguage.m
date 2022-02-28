function true = isRxnWrittenInLanguage(model, posRxn, languageToCheck, database)

true = 0;
posMets = find(model.S(:,posRxn));
mets = model.mets(posMets);
if allMetsInlanguage(mets,languageToCheck, database)
    true = 1;
end

end