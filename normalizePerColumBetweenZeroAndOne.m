function userDataset = normalizePerColumBetweenZeroAndOne(userDataset)


a = userDataset;
x = [min(a,[],1);max(a,[],1)];
b = bsxfun(@minus,a,x(1,:));
userDataset = bsxfun(@rdivide,b,diff(x,1,1));

end