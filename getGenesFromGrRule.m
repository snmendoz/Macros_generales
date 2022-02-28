function genes = getGenesFromGrRule(grRule)

genes = unique(splitString(regexprep(grRule,{'\(|\)','\ or\ |\ and\ |\ &\ |\ \|\ '},{'',' '})));

end