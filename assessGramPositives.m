function assessGramPositives
cd('D:\Dropbox\Databases\BIGG')
models = {'iAF692';'iHN637';'iNF517';'iNJ661';'iSB619';'iYO844'};

for i =1:length(models)
   load(models{i})
   model = eval(models{i});
   c = model.mets(find(cellfun(@isempty, strfind(model.mets,'_c'))==0));
   p = model.mets(find(cellfun(@isempty, strfind(model.mets,'_p'))==0));
   
   if isempty(c)
       disp('')
   else
       if ~isempty(p)
           disp('')
       end
   end
end

end