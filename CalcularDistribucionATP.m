function [h,s,ATP_het,ATP_syn]=CalcularDistribucionATP(model,solution)
[~,a,~]=intersect(model.rxns,{'ACETATEKIN_RXN','PEPDEPHOS_RXN','PHOSGLYPHOS_RXN'});
[~,a4,~]=intersect(model.rxns,{'GLUCOKIN_RXN','FRUCTOKINASE_RXN','TRANS_RXNK9E_535','TRANS_RXNK9E_536'});
[~,a3,~]=intersect(model.rxns,{'ATPSYN_RXN'});
ATP_het=sum(abs(solution(a)))-sum(abs(solution(a4)));
ATP_syn=abs(solution(a3));

total=ATP_het+ATP_syn;

h=ATP_het/total;
s=ATP_syn/total;

end