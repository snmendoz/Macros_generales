function [model, mets ] = homogenizeModelAndMets(model, mets, desiredFormat)

switch desiredFormat
    case 'cobra'
        if allMetsInCBMPYFormat(model.mets) 
            model = transformModelToCOBRAFormat(model);
        end
        if allMetsInCBMPYFormat(mets) 
            mets = transformMetsToCOBRAFormat(mets);
        end 
    case 'cbmpy'
        if allMetsInCOBRAFormat(model.mets) 
            model = transformModelToCBMPYFormat(model);
        end
        if allMetsInCOBRAFormat(mets) 
            mets = transformMetsToCBMPYFormat(mets);
        end 

end