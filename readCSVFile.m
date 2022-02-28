function data = readCSVFile(filename, separator)

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
    data{n,1} = line{1}; 
    data{n,2} = line{2};
    tline = fgetl(fid);
end

data = data(1:n,:);

end