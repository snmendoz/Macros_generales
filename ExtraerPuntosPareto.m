function [x,y,x_adicionales,y_adicionales]=ExtraerPuntosPareto(Dir,File)

current=pwd;
cd(Dir);

fi = fopen (File,'r+');
x=[];
y=[];
tline = fgetl(fi);
while ~isempty(tline)
    numbers=regexp(tline,' ','split');
    if strcmp(numbers{1},'0')
    elseif strcmp(numbers{1},'1')
        x=[x;str2num(numbers{2})];
        y=[y;str2num(numbers{3})];
    end
    tline = fgetl(fi);
    if tline==-1
        break;
    end
end
fclose(fi);

[x_ord,ind]=sort(x);
y_ord=y(ind);
x_adicionales=[];
y_adicionales=[];
for i=1:length(x_ord)-1
    x_adicionales=[x_adicionales;mean([x_ord(i),x_ord(i+1)])];
    
    x1=x_ord(i);
    x2=x_ord(i+1);
    y1=y_ord(i);
    y2=y_ord(i+1);
    x_i=mean([x_ord(i),x_ord(i+1)]);
    
    m=(y2-y1)/(x2-x1);
    y_i=m*(x_i-x1)+y1;
    y_adicionales=[y_adicionales;y_i];
end
x=[x_ord;x_adicionales];
[x_ord_b,ind]=sort(x);
y=[y_ord;y_adicionales];
y_ord_b=y(ind);

x=x_ord;
y=y_ord;
x_adicionales=x_ord_b;
y_adicionales=y_ord_b;
cd(current);

end