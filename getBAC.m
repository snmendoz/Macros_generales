function BAC = getBAC(true, pred)

TP = length(intersect(find(true==1),find(pred==1)));
TN = length(intersect(find(true==2),find(pred==2)));
FP = length(intersect(find(true==2),find(pred==1)));
FN = length(intersect(find(true==1),find(pred==2)));
SEN = TP/(TP + FN);
SPE = TN/(TN + FP);
BAC = 0.5*(SEN+SPE);

end