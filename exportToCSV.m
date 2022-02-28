function exportToCSV(fileName, matrix)

fid = fopen([fileName '.csv'],'w+');
for i = 1:size(matrix,1)
    if any(~cellfun(@(x) ischar(x), matrix(i,:)))
        for j = 1:length(matrix(i,:))
            if ~ischar(matrix{i,j})
                matrix{i,j} = num2str(matrix{i,j});
            end
        end
    end
    fprintf(fid, [strjoin(matrix(i,:),',') '\n']);      
end

fclose(fid);


end