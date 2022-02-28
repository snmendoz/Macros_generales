function checkAuReMeTemplates(path)

load('D:\Dropbox\Databases\BIGG\bigg_85.mat')

cd(path);
if exist('rxnRep.xls','file') ==2
    delete('rxnRep.xls')
end
if exist('metRep.xls','file') ==2
    delete('metRep.xls')
end


data = dir();
folders = {};
for i = 1:length(data)
    if ~strcmp(data(i).name,'.') && ~strcmp(data(i).name,'..')
        folders = [folders; data(i).name];
    end
end

rxns = cell(size(folders));
mets = cell(size(folders));
for i = 1:length(folders)
    cd([path filesep folders{i}])
    model_i = readCbModel('metabolic_model.sbml');
    model_i.mets = regexprep(model_i.mets,'\[(.)\]','_$1');
    
    rxnDiff = setdiff(model_i.rxns,bigg.rxns);
    metDiff = setdiff(model_i.mets,bigg.mets);
    if ~isempty(rxnDiff)
        rxns{i} = rxnDiff;
    end
    if ~isempty(rxnDiff)
        mets{i} = metDiff;
    end
end
cd(path);
for i = 1:length(rxns)
    if ~isempty(rxns{i})
        xlswrite('rxnRep', rxns{i},folders{i});
    end
    if ~isempty(mets{i})
        xlswrite('metRep', mets{i},folders{i});
    end
end
    
end