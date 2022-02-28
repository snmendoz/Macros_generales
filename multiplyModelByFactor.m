function modelOut = multiplyModelByFactor(model, factor)

modelOut = model;
modelOut.S = model.S*factor;
modelOut.lb = model.lb.*factor;
modelOut.ub = model.ub.*factor;
modelOut.b = model.b.*factor;
end