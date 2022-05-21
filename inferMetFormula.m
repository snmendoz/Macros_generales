function [model, elements, realCoefs] = inferMetFormula(model,rxn, aprox)

if nargin <3
   aprox = 0; 
end
realCoefs = [];

mets = getMetaboliteIDsFromRxns(model, {rxn});
mets = mets{1};
missingFormula = find(cellfun(@isempty, (model.metFormulas(getPosOfElementsInArray(mets, model.mets)))));

posMissingMet = getPosOfElementsInArray(mets(missingFormula), model.mets);
posRxn = getPosOfElementsInArray({rxn}, model.rxns);

elements = {'C','H','O','N','P','S','Na','Mg','Cl','K','Ca','Mn','Fe','Ni','Co','Cu','Zn','As','Se','Ag','Cd','W','Hg','Mo','I','R','X'};
numbers = zeros(size(elements));

rest = setdiff(mets, mets(missingFormula));
coeff_missing_met = full(model.S(posMissingMet, posRxn));

for k = 1:length(elements)
    disp(k)
    new_coeff = 0;
    for i = 1:length(rest)
        disp(i)
        pos_rest_i = getPosOfElementsInArray(rest(i), model.mets);
        coeff_i = full(model.S(pos_rest_i,posRxn));
        formula_i = model.metFormulas{pos_rest_i};
        n = numAtomsOfElementInFormula(formula_i, elements{k});
        
        if coeff_i>0 && coeff_missing_met>0 || coeff_i<0 && coeff_missing_met<0
             new_coeff = new_coeff - abs(coeff_i)*n;
        else
            new_coeff = new_coeff + abs(coeff_i)*n;
        end
        
    end
    new_coeff = new_coeff/abs(coeff_missing_met);
    realCoefs = [realCoefs;new_coeff];
    if aprox
        new_coeff = round(new_coeff);
    end
    numbers(k)=new_coeff;
end

notZero = find(numbers~=0);
new_formula = createFormulaFromElements(elements(notZero), numbers(notZero));

model.metFormulas{posMissingMet} = new_formula;

new_charge = 0;
for i = 1:length(rest)
    pos_rest_i = getPosOfElementsInArray(rest(i), model.mets);
    charge_i = model.metCharges(pos_rest_i);
    coeff_i = full(model.S(pos_rest_i,posRxn));
    
    if coeff_i>0 && coeff_missing_met>0 || coeff_i<0 && coeff_missing_met<0
        new_charge = new_charge - charge_i*abs(coeff_i);
    else
        new_charge = new_charge + charge_i*abs(coeff_i);
    end
    
end
format short;
new_charge = new_charge/abs(coeff_missing_met);
if aprox
        new_charge = round(new_charge);
end
model.metCharges(posMissingMet) = new_charge;

end