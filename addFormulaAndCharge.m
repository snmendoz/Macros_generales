function addFormulaAndCharge(database)
cd('D:\Dropbox\Databases\BIGG')
load(database)
formulas = cell(size(bigg.mets));
charges = nan(size(bigg.mets));
posRevisar = zeros(1000,1);
posDoble = zeros(1000,1);
posFormulaVacia = zeros(1000,1);
posChargeVacia = zeros(1000,1);

cont_posRevisar = 0;
cont_posDoble = 0;
cont_posFormulaVacia = 0;
cont_posChargeVacia = 0;

mets = regexprep(bigg.mets,'(.*)\[.*\]$' ,'$1');
mets = unique(mets);

for i = 1:length(mets)
    disp(i)
    pos = find(strncmp(bigg.mets, [mets{i} '['], length([mets{i} '['])));
    [s, r] = system(['curl http://bigg.ucsd.edu/api/v2/universal/metabolites/' mets{i}]);
    if s==0
        posF = strfind(r, 'formulae');
        if ~isempty(posF)
            posComma = strfind(r, ',');
            posComma1 = posComma(posComma>posF);
            r1 = r(posF:posComma1(1));
            posSqrBr1 = strfind(r1, '[');
            posSqrBr2 = strfind(r1, ']');
            if isempty(posSqrBr2)
                formula = regexprep(r1(posSqrBr1(1)+1: end-2), '"', '');
            else
                formula = regexprep(r1(posSqrBr1(1)+1: posSqrBr2(1)-1), '"', '');
            end
            
            if isempty(formula)
                cont_posFormulaVacia = cont_posFormulaVacia + 1;
                posFormulaVacia(cont_posFormulaVacia) = i;
            else
                if isValidFormula(formula)
                    formulas(pos) = {formula};
                else
                    disp('')
                end
            end
        end
        
        posC = strfind(r, 'charges');
        if ~isempty(posC)
            r2 = r(posC:end);
            posSqrBr1 = strfind(r2, '[');
            posSqrBr2 = strfind(r2, ']');
            charge = r2(posSqrBr1(1)+1: posSqrBr2(1)-1);
            
            if isempty(charge)
                cont_posChargeVacia = cont_posChargeVacia + 1;
                posChargeVacia(cont_posChargeVacia) = i;
            else
                if strfind(charge, ',')
                    cont_posDoble = cont_posDoble + 1;
                    posDoble(cont_posDoble) = i;
                else
                    charges(pos) = str2double(charge);
                end
            end         
        end                     
    else
        disp('')
        cont_posRevisar = cont_posRevisar + 1;
        posRevisar(cont_posRevisar) = i;
    end
end

posDoble = posDoble(1:cont_posDoble);
posRevisar = posRevisar(1:cont_posRevisar);
posChargeVacia = posChargeVacia(1:cont_posChargeVacia);
posFormulaVacia = posFormulaVacia(1:cont_posFormulaVacia);

save('posDoble', 'posDoble')
save('posChargeVacia', 'posChargeVacia')
save('posFormulaVacia', 'posFormulaVacia')
save('posRevisar', 'posRevisar')

bigg.metFormulas = formulas;
bigg.metCharges = charges;
save(database, 'bigg');

end