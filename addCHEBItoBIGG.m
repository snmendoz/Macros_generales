function addCHEBItoBIGG(database)
cd('D:\Dropbox\Databases\BIGG')
load(database)
CHEBI = cell(size(bigg.mets));
posRevisar = zeros(1000,1);
cont_posRevisar = 0;
mets = regexprep(bigg.mets,'(.*)_(.*)$' ,'$1');
for i = 1:length(bigg.mets)
    disp(i)
    [s, r] = system(['curl http://bigg.ucsd.edu/api/v2/universal/metabolites/' mets{i}]);
    if s==0
        posF = strfind(r, '"CHEBI:');
        chebis = {};
        if ~isempty(posF)
            for j = 1:length(posF)
                posComma = strfind(r, '"');
                posComma1 = posComma(posComma>posF(j));
                r1 = r(posF(j):posComma1(1));
                chebis = [chebis, regexprep(r1,{'"','CHEBI:'},{'', ''})];
            end
        end
        CHEBI{i} = strjoin(chebis,','); 
        
    else
        disp('')
        cont_posRevisar = cont_posRevisar + 1;
        posRevisar(cont_posRevisar) = i;
    end
end

posRevisar = posRevisar(1:cont_posRevisar);
save('posRevisar', 'posRevisar')

bigg.metChEBIID = CHEBI;
for i = 1:length(bigg.metChEBIID); if isempty(bigg.metChEBIID{i}); bigg.metChEBIID{i} =''; end; end;
save(database, 'bigg');

end