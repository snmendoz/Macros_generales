function massBalances = getMassEquationsFromModel_from_met(model, pos)

if ischar(pos)
   pos = getPosOfElementsInArray({pos}, model.mets); 
end

massBalances = cell(length(pos),1);
variables = strcat('v_', model.rxns);
for i = 1:length(pos)
     pos_variables_equation_i = find(model.S(pos(i),:));
     coefs = full(model.S(pos(i),pos_variables_equation_i));
     equation_i = [model.mets{pos(i)} ': '];
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
     if strcmp('E', model.csense(pos(i)))
         equation_i = [equation_i ' = ' num2str(model.b(pos(i)))];
     elseif strcmp('L', model.csense(pos(i)))
         equation_i = [equation_i ' <= ' num2str(model.b(pos(i)))];
     else
         equation_i = [equation_i ' >= ' num2str(model.b(pos(i)))];
     end
     massBalances{i} = equation_i;
end


end