function plotVennForThreeModels(models, opt)
figure;

if opt==1
    s1 = models{1}.rxns;
    s2 = models{2}.rxns;
    s3 = models{3}.rxns;
else
    s1 = models{1}.mets;
    s2 = models{2}.mets;
    s3 = models{3}.mets;
end

%uvarum
% figure
int = intersect(s1, s2);
% A = [length(s1) length(s2)];
% [H, S] = venn(A,length(int),'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');

% figure
% int = intersect(s1, s3);
% A = [length(s1) length(s3)];
% [H, S] = venn(A,length(int),'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');

int2 = intersect(int, s3);

% [i12 i13 i23 i123]
A = [length(s1) length(s2) length(s3)];
I = [length(intersect(s1, s2)),...
    length(intersect(s1, s3)),...
    length(intersect(s2, s3))...
    length(int2)];

try
    [H, S] = venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black');
catch
    warning('issue with set function')
end

end