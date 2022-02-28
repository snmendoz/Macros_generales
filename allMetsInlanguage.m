function true = allMetsInlanguage(mets,language, database)

true = 0;
if all(cellfun(@(x) isMetInLanguage(x,language,database), mets))
    true = 1;
end

end