function relaxedModel = relaxModel(model, relaxedMets)

relaxedModel = model; 
for i = 1:length(relaxedMets)
   relaxedModel = addReaction(relaxedModel, [relaxedMets(i) '_DM'], relaxedMets(i), -1, 0, 0, 1000);  
end

end