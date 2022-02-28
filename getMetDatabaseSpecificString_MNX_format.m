function outputMet = getMetDatabaseSpecificString_MNX_format(met,language)

met_woc = removeCompartmentFromMets({met});
outputMet = [language, ':', met_woc{1}];

end