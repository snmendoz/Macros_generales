function coef = findUptakeRateUsingFitting(xs,ys, mu, ...
    initialCompoundConcentration, initialBiomass, type)

%we assume this formula C(t) = C_0 + (rate/mu)*X_0*(exp(mu*t) - 1).
%we want to estimate rate

opt = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',0,...
    'Upper',100,...
    'Startpoint',0.01);
if strcmp(type,'P')
    f = fittype([num2str(initialCompoundConcentration) ' + (a/' num2str(mu) ')*' num2str(initialBiomass) '*(exp(' num2str(mu) '*x)-1)'],'options',opt);
elseif strcmp(type,'N')
    f = fittype([num2str(initialCompoundConcentration) ' - (a/' num2str(mu) ')*' num2str(initialBiomass) '*(exp(' num2str(mu) '*x)-1)'],'options',opt);
end
[coef,~] = fit(xs,ys,f);
% figure
% plot(xs,ys,'bo')
% hold on
% plot(coef,'m')

end