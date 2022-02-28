function newFormula = createFormulaFromElements(elements, numbers)

newFormula = '';
format longG
for i = 1:length(elements)
    newFormula = [newFormula elements{i} num2str(abs(numbers(i)))];
    
end

end