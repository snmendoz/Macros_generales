function cobra2vlp(model, objInd, outputFileName)

vlp = cobra2prob(model, objInd);
prob2vlp(vlp, [outputFileName '.vlp']);

end