function [threshold,decimalPecision] = findLimitingNutrientThreshold(model, nutrient,seedLB, seedUB)
maxAbsoluteIterations = 10;
maxRelativeIterations = 6;
if nargin < 3
    seedLB = 1;
end
if nargin < 4
    seedUB = 10;
end
currentStep = 1;
currentRelativeStep = 0;

is = IsNutrientLimiting(model, nutrient);

if is
    while currentStep<maxAbsoluteIterations && currentRelativeStep<maxRelativeIterations 
        array = linspace(seedLB,seedUB,10); 
        mus = zeros(size(array));
        for i = 1:length(array)
            modelAux = model;
            modelAux = changeRxnBounds(modelAux,nutrient, -array(i),'l');
            fba = optimizeCbModel(modelAux);
            mus(i) = fba.f;
        end
        if currentRelativeStep>0
            if floor(mus(1)*10^7)/10^7 == floor(mus(end)*10^7)/10^7
                break;
            end
        end
        mus = floor(mus*10^8)/10^8;
        
        if mus(1) == mus(end)
            seedLB = array(1)/10;
            seedUB = array(end)/10;
        elseif length(find(mus==mus(end)))>1
            currentRelativeStep = currentRelativeStep+1;
            pos = find(mus==mus(end));
            seedLB = array(pos(1)-1);
            seedUB = array(pos(1)); 
        elseif mus(end)>mus(end-1)
            currentRelativeStep = currentRelativeStep+1;
            seedLB = array(end-1);
            seedUB = array(end);
        end
        
        currentStep = currentStep+1;
    end
end

posRound = find(arrayfun(@(x) floor(seedLB*10^x)/10^x~=floor(seedUB*10^x)/10^x, 1:14));
posRound = posRound(1);
threshold = floor(seedLB*10^(posRound-1))/10^(posRound-1);
decimalPecision = posRound-1;
end