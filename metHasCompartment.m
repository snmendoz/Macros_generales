function has = metHasCompartment(met, listOfPossibleCompartments)

has = any(cellfun(@(x) ~isempty(regexp(met, ['_' x '$|[' x ']$'], 'once')), listOfPossibleCompartments));
end