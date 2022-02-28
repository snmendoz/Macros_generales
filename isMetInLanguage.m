function true = isMetInLanguage(met,language,database)

true = 0;
if ismember(getMetDatabaseSpecificString_MNX_format(met,language),database)
    true = 1;
end

end