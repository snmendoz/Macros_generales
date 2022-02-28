function [rxns, pos] = getReactionsInvolvingMetaboliteIdentifiers(model, identifiers)
rxns = {};
pos = [];
positions = cell(length(identifiers),1);
for i = 1:length(identifiers)
    possibleMets = getMetabolitesWithCompartmentsInModelFromIdentifier(model, identifiers{i});
    [~, ~, ~, ~, positions_i] = getRxnsFromMets(model, possibleMets);
    positions{i} = positions_i;
end
final_positions = positions{1};
for i = 2:length(identifiers)
    final_positions = intersect(final_positions, positions{i});
end
if ~isempty(final_positions)
    rxns =model.rxns(final_positions);
    pos = final_positions;
end

end