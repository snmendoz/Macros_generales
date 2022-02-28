function Info = compareStoic(model1,model2,rxnIDs)

if nargin<3
    [~,ind1,ind2]=intersect(model1.rxns,model2.rxns);
else
    for i=1:length(rxnIDs)
    [~,pos1]=ismember(model1.rxns,rxnIDs(i));
    [~,pos2]=ismember(model2.rxns,rxnIDs(i));
    ind1(i)=find(pos1);
    ind2(i)=find(pos2);
    end
end

same_mets=zeros(length(ind1),1);
same_stoic=zeros(length(ind1),1);
same_stoic_rev=zeros(length(ind1),1);
same_lb_rev=zeros(length(ind1),1);
same_ub_rev=zeros(length(ind1),1);
bloq1=model1.lb==0 & model1.ub==0;
forward1=model1.lb>=0 & model1.ub>=0 & ~bloq1;
backward1=model1.lb<=0 & model1.ub<=0 & ~bloq1;
rev1=model1.lb<0 & model1.ub>0; 
bloq2=model2.lb==0 & model2.ub==0;
forward2=model2.lb>=0 & model2.ub>=0 & ~bloq2;
backward2=model2.lb<=0 & model2.ub<=0 & ~bloq2;
rev2=model2.lb<0 & model2.ub>0; 

same_rev_overall=ones(length(ind1),1);

for i=1:length(ind1)
    pos_mets1=find(model1.S(:,ind1(i)));
    pos_mets2=find(model2.S(:,ind2(i)));
    mets1=model1.mets(pos_mets1);
    mets2=model2.mets(pos_mets2);
    
    if length(mets1)==length(mets2) && length(intersect(mets1,mets2))==length(mets1)
        same_mets(i)=1;
        stoic_i=1;
        stoic_rev_i=1;
        for j=1:length(mets1)
            pos1=pos_mets1(j);
            pos2=strmatch(mets1(j),model2.mets,'exact');
            stoic1=full(model1.S(pos1,ind1(i)));
            stoic2=full(model2.S(pos2,ind2(i)));
            if stoic1~=stoic2; stoic_i=0; end; 
            if stoic1~=-stoic2; stoic_rev_i=0; end; 
        end
        if stoic_i;
            same_stoic(i)=1;
            if model1.lb(ind1(i))==model2.lb(ind2(i))
                same_lb_rev(i)=1;
            end
            if model1.ub(ind1(i))==model2.ub(ind2(i))
               same_ub_rev(i)=1;
            end 
        else
            same_stoic(i)=0;
        end
        if stoic_rev_i;
            same_stoic_rev(i)=1;
            if model1.lb(ind1(i))==-model2.ub(ind2(i))
                same_lb_rev(i)=1;
            end
            if model1.ub(ind1(i))==-model2.ub(ind2(i))
                same_ub_rev(i)=1;
            end
            if (rev2(ind2(i)) && ~rev1(ind1(i))) || (~rev2(ind2(i)) && rev1(ind1(i)))
                same_rev_overall(i)=0;
            end
            
        else
            same_stoic_rev(i)=0;
        end
    else
        same_mets(i)=0;
        same_stoic(i)=0;
        %metabolites do not match
    end  
end
same_stoic_overall=same_stoic | same_stoic_rev;
n_same_stoic_overall=length(find(same_stoic_overall));
different_stoic_overall=ind1(find(same_stoic_overall==0));

Info=[];
for i=1:length(different_stoic_overall)
    eq1=getRxn_cobraFormat(model1,ind1(find(ind1==different_stoic_overall(i))));
    eq2=getRxn_cobraFormat(model2,ind2(find(ind1==different_stoic_overall(i))));
    Info=[Info; [[{'rxnID'};{'model1'};{'model2'}],[model1.rxns(different_stoic_overall(i));eq1;eq2]]];
end

eq1=getRxn_cobraFormat(model1,ind1(find(same_ub_rev==0)));
eq2=getRxn_cobraFormat(model2,ind2(find(same_ub_rev==0)));

lb1=model1.lb(ind1(find(same_ub_rev==0)));
ub1=model1.ub(ind1(find(same_ub_rev==0)));
lb2=model2.lb(ind2(find(same_ub_rev==0)));
ub2=model2.ub(ind2(find(same_ub_rev==0)));

Info2=[eq1,eq2,num2cell(lb1),num2cell(ub1),num2cell(lb2),num2cell(ub2)];
end