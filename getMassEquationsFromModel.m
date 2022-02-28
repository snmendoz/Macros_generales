function massBalances = getMassEquationsFromModel(model,fileName)

massBalances = cell(length(model.mets),1);
variables = strcat('v_', model.rxns);
for i = 1:length(model.mets)
     pos_variables_equation_i = find(model.S(i,:));
     coefs = full(model.S(i,pos_variables_equation_i));
     equation_i = [model.mets{i} ': '];
     for j = 1:length(pos_variables_equation_i)
         if j ==1
             equation_i = [equation_i num2str(coefs(j)) ' ' variables{pos_variables_equation_i(j)}];
         else
             if coefs(j)>0
                 equation_i = [equation_i ' + ' num2str(coefs(j)) ' ' variables{pos_variables_equation_i(j)}];
             else
                 equation_i = [equation_i ' - ' (num2str(abs(coefs(j)))) ' ' variables{pos_variables_equation_i(j)}];
             end
         end
     end
     if strcmp('E', model.csense(i))
         equation_i = [equation_i ' = ' num2str(model.b(i))];
     elseif strcmp('L', model.csense(i))
         equation_i = [equation_i ' <= ' num2str(model.b(i))];
     else
         equation_i = [equation_i ' >= ' num2str(model.b(i))];
     end
     massBalances{i} = equation_i;
end

if nargin>1
    xlswrite(fileName, massBalances)
end

end