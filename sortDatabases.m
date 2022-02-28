function [trainingDataset, new_testDataset, features] = sortDatabases(trainingDataset,testDataset,traning_features, test_features)

features = traning_features;
new_testDataset = zeros(size(testDataset,1),size(trainingDataset,2));

for i = 1:length(traning_features)
    pos = getPosOfElementsInArray(traning_features(i), test_features);
    if ~isempty(pos)
        new_testDataset(:,i) = testDataset(:,pos);
    end
end


end