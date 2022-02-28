function checkBigg(database)
load(database);
Redundantes=zeros(length(bigg.rxns));
posRedundantes=cell(length(bigg.rxns),1);

for i=10:length(bigg.rxns)
    [Redundant,posRedundants]=IsDuplicate(bigg,i);
    if Redundant
        Redundantes(i)=1;
        posRedundantes{i}=[i;posRedundants];
    end
end
Redundantes_pos=find(Redundantes);
posRedundantes=posRedundantes(Redundantes_pos);

end