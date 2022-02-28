function data = readCSVFile_general(filename, separator)

if nargin<2
   separator = ','; 
end

fid = fopen(filename,'r+');
tline = fgetl(fid);
data = cell(5000,2);
n = 0;
while ischar(tline)
    line = strsplit(tline, separator);
    n = n+1;
    for j = 1:length(line)
        data{n,j} = line{j};
    end
    tline = fgetl(fid);
end

fclose(fid);
data = data(1:n,:);

end