function coef = fitPh(xs, ys, weights)

if nargin<3
    weights = ones(size(xs));
end

opt = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,0, 0, 0],...
    'Upper',[Inf,Inf,Inf, Inf],...
    'Startpoint',[1 1 1,1],...
    'Robust','on',...
    'Weights',weights);
f = fittype('-(X0*K/((K-X0)*exp(-mu*x)+X0))+beta','options',opt);

[coef,~] = fit(xs,ys,f);


end