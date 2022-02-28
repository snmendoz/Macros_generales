function [type,x,Adj]=read_ParetoFront(ImgfileName,AdjFileName,Directory)

if nargin>2 && ~isempty(Directory)
    current=pwd;
    cd(Directory)
end

fi = fopen (ImgfileName,'r+');
x=[];
type=[];
tline = fgetl(fi);
while ~isempty(tline)
    numbers=regexp(tline,' ','split');
    if strcmp(numbers{1},'0')
        type=[type;0];
        for i=2:length(numbers)
            num(i-1)=str2num(numbers{i});
            
        end
        x=[x;num];
    elseif strcmp(numbers{1},'1')
        type=[type;1];
        for i=2:length(numbers)
            num(i-1)=str2num(numbers{i});
            
        end
        x=[x;num];
    end
    tline = fgetl(fi);
    if tline==-1
        break;
    end
end
fclose(fi);


fi2 = fopen (AdjFileName,'r+');
tline = fgetl(fi2);
Adj=[];

while ~isempty(tline)
    numbers=regexp(tline,' ','split');
    num=[];
    for i=1:length(numbers)
        num(i)=str2num(numbers{i});
    end
    Adj=[Adj;{num}];
    tline = fgetl(fi2);
    if tline==-1
        break;
    end
end
fclose(fi2);

if nargin>2 && ~isempty(Directory)
    cd(current);
end
end