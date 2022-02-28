function WriteReactionsToAddInAuremeFormat(fileName, forAdding, model)

%agregadas desde bigg
fi=fopen([ fileName '.txt'],'w+');
for i = 1:length(forAdding)
    pos = find(strcmp(model.rxns,forAdding(i)));
    fprintf(fi,[ 'reaction_id\t' forAdding{i} '\n']);
    fprintf(fi,'comment\tmanual curation\n');
    if model.lb(pos)<0
        fprintf(fi,'reversible\ttrue\n');
    else
        fprintf(fi,'reversible\tfalse\n');
    end
    fprintf(fi,[ 'linked_gene\t' model.grRules{pos} '\n']);
    
    pos_mets_rxn_i=find(model.S(:,pos));
    coef = full(model.S(pos_mets_rxn_i,pos));
    [~, ind] = sort(coef);
    pos_mets_rxn_i = pos_mets_rxn_i(ind);
    for k=1:length(pos_mets_rxn_i)
        comp = model.mets{pos_mets_rxn_i(k)};
        comp = comp(end-1);
        met=regexprep(model.mets{pos_mets_rxn_i(k)},{'\[c\]','\[e\]','\[e\]'},{'','',''});
        if full(model.S(pos_mets_rxn_i(k),pos))<0
            fprintf(fi,[ 'reactant\t%3.1f:' met ':' comp '\n'],abs(full(model.S(pos_mets_rxn_i(k),pos))));
        else
            fprintf(fi,[ 'product\t%3.1f:' met ':' comp '\n'],abs(full(model.S(pos_mets_rxn_i(k),pos))));
        end
    end
    fprintf(fi,'\n');    
end
fclose(fi);

end