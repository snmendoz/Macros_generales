function exportInputsForEMAF(modelIrrev)

posEX = find(~cellfun(@isempty, strfind(modelIrrev.rxns, 'EXi')));
rxns = modelIrrev.rxns(posEX);
data = [{'Reaction ID'}; rxns];
fid = fopen('all.csv','wt');
if fid>0
    for k=1:size(data,1)
        fprintf(fid,'%s\n',data{k,:});
    end
    fclose(fid);
end

% csvwrite('all.csv',) 

changeCobraSolver('gurobi')
biomass = modelIrrev.rxns(find(modelIrrev.c));
fba = optimizeCbModel(modelIrrev);
lb = fba.f*0.6;
ub = fba.f;

data = [[{'Reaction ID'},{'LB'},{'UB'}];biomass, num2str(lb,'%1.2f'), num2str(ub,'%1.2f')];

fid = fopen('basic.csv','wt');
if fid>0
    for k=1:size(data,1)
        fprintf(fid,'%s,%s,%s\n',data{k,:});
    end
    fclose(fid);
end

% csvwrite('basic.csv',[[{'Reaction ID'},{'LB'},{'UB'}];biomass, lb, ub]) 

end