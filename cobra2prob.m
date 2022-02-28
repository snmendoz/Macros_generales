function [vlp,P,B,a,b,l,s,opt_dir,Y,Z,c]=cobra2prob(cobra,obj)
% We are using bensolve-2.0.1:
% B is coefficient matrix
% P is objective Marix
% a is lower bounds for B
% b is upper bounds for B
% l is lower bounds of variables
% s is upper bounds of variables
% opt_dir is direction: 1 min, -1 max
% Y,Z and c are part of cone definition. If empty => MOLP
    [m, n] = size(cobra.S);
    q = length(obj);
    vlp.B = cobra.S;
    vlp.a = zeros(m,1);
    vlp.b = zeros(m,1);
    vlp.l = cobra.lb;
    vlp.s = cobra.ub;
    vlp.P = zeros(q,n);
    vlp.opt_dir=-1;
    for i=1:q
        vlp.P(i,obj(i))=1;
    end
    vlp.Y=[]; vlp.Z=[];vlp.c=[];
    
    B=vlp.B ;
    a=vlp.a ;
    b=vlp.b ;
    l=vlp.l ;
    s=vlp.s ;
    P=vlp.P ;
    opt_dir=vlp.opt_dir;
    Y=[]; Z=[];c=[];
    
    
end