function writeCSVFile_general(array, fileName,separator)

if nargin<3
   separator = ','; 
end

if ~isempty(strfind(fileName,'.txt'))
    fileID = fopen(fileName,'wt');
else
    fileID = fopen([fileName '.txt'],'wt');
end

if isnumeric(array)
    array2 = cell(size(array));
    for i = 1:size(array,1)
        for j = 1:size(array,2)
            array2{i,j} = num2str(array(i,j));
        end
    end
else
    array2 = array;
end
array = array2;

for i=1:size(array,1)
    fprintf(fileID, '%s\n', strjoin(array(i,:),separator));
end
fclose(fileID);


end