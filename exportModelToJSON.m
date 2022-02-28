function exportModelToJSON(inputName, outputName, folder)

exportModelToJSON_filepath = which('exportModelToJSON.py');
system(['python "' exportModelToJSON_filepath '" "' inputName '" "' outputName '" "' folder '"'])

end

