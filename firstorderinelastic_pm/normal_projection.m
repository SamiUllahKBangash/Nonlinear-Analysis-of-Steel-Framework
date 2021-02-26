function [msol,psol,alphasol]=normal_projection(m1,p1,Py,Mp)

%Welcome to drift corrector routine. This script is written to return trial
%fint values (that have been scaled to remain within expanded surface) back
%to original yield surface. Corrected fint trial values are returned via
%gradient to the original surface as opposed to radial return used
%previously. Matlab's vpasolve is invoked to solve three highly nonlinear
%equations for m, p, alpha parameters. 

      
syms m p alp real
eq1=m*(2*alp + 7*alp*p^2 + Mp)== m1*Mp;

eq2=p*(2*alp + 7*alp*m^2 + Py) == p1*Py;

eq3=m^2 + p^2 + 3.5*p^2*m^2==1;

%[sol1,sol2,sol3]=vpasolve(eq1,eq2,eq3,[m,p,alp]);
%S=vpasolve(eq1,eq2,eq3,m,[0 1],p,[0 1]);
vars=[m p alp];

range=[0 sign(m1)*1;0 sign(p1)*1;NaN NaN];
S=vpasolve(eq1,eq2,eq3,vars,range);
msol=S.m;
psol=S.p;
alphasol=S.alp;
%{
I=find((S.m)*m1>0 & (S.p)*p1>0 );
mcorrec=S.m(I);
pcorrec=S.p(I);
%alpha=S.alp(I)
%}
end
