function getRxnBalance(model,pos)
info = cell(100,6);
%substrate
[left, pos_left] = getMetabolitesInLeftFromRxn(model,pos);
[right, pos_rigth] = getMetabolitesInRigthFromRxn(model,pos);


n_left = length(left);
n_right = length(right);
for i =1:n_left
    info{i,1} = left{i};
    info{i,2} = model.metFormulas{pos_left(i)};
    info{i,3} = model.metCharges(pos_left(i));
end

for i =1:n_right
    info{i,4} = right{i};
    info{i,5} = model.metFormulas{pos_rigth(i)};
    info{i,6} = model.metCharges(pos_rigth(i));
end



end