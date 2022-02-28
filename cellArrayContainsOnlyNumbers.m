function isTrue = cellArrayContainsOnlyNumbers(array)

isTrue = all(cellfun(@(x) ~isnan(str2double(x)), array));

end