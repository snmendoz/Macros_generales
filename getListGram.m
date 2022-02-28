function [bacteriaList, gramPosList, gramNegList, bacteria, gramPos, gramNeg] = getListGram()
load('D:\Projects\BIGG\models_refined.mat')
bact_bool = ones(size(models_refined));
bact_grampos_bool = zeros(size(models_refined));
bact_gramneg_bool = zeros(size(models_refined));

bacteriaList = {};
gramPosList = {};
gramNegList = {};
for i = 1:length(models_refined)
    model_i = models_refined{i};
    [tok,rem] = strtok(model_i.mets,'[');
    compartments = unique(regexprep(rem,{'[',']'},{'',''}));

    if ~isempty(setdiff(compartments,{'c','p','e'}))
        bact_bool(i)=0;
    end
    if bact_bool(i)==1; bacteriaList = [bacteriaList; model_i.description]; end
    if length(intersect({'c','e'}, compartments)) == length(compartments) && length(compartments) == 2
        bact_grampos_bool(i)=1;
        gramPosList = [gramPosList; model_i.description];
    end
    if length(intersect({'c','e','p'}, compartments)) == length(compartments) && length(compartments) == 3
        bact_gramneg_bool(i)=1;
        gramNegList = [gramNegList; model_i.description];
    end
end

bacteria = models_refined(find(bact_bool));
gramPos = models_refined(find(bact_grampos_bool));
gramNeg = models_refined(find(bact_gramneg_bool));

end