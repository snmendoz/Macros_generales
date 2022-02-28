function charge = calculateChargeFromReference(userFormula, refFormula, refCharge)

userFormula = strrep(userFormula, 'FULLR3', 'R3');
userFormula = strrep(userFormula, 'FULLR2', 'R2');
userFormula = strrep(userFormula, 'FULLR', 'R');

n_userFormula = numAtomsOfElementInFormula(userFormula, 'H');
n_refFormula = numAtomsOfElementInFormula(refFormula, 'H');

charge = refCharge + (n_userFormula-n_refFormula);

end