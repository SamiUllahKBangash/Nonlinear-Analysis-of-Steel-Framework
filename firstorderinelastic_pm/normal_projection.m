function [msol,psol,alphasol]=normal_projection(m1,p1,Py,Mp)

%Welcome to drift corrector routine. This script is written to return trial
%fint values (that have been scaled to remain within expanded surface) back
%to original yield surface. Corrected fint trial values are returned via
%gradient to the original surface as opposed to radial return used
%previously. Matlab's vpasolve is invoked to solve three highly nonlinear
%equations for m, p, alpha parameters. 



%{
m1=0.4;
p1=0.3;
Py=955; %kips/K
Mp=4840; %Kin
mi=sign(m1)*0.5;
converge=false;
j=1;
while ~converge && j<=100
    
    if p1>0
        m=((2+7*mi^2)*Mp*(m1 + 3.5*m1*mi^2 - mi + 3.5*mi^3))/(9*Py*(p1*sqrt(1+3.5*mi^2) - sqrt(1-mi^2)));
    else
        m=((2+7*mi^2)*Mp*(m1 + 3.5*m1*mi^2 - mi + 3.5*mi^3))/(9*Py*(p1*sqrt(1+3.5*mi^2) - sqrt(1-mi^2)));
    end
    
    if (abs(m-mi)/mi)*100 <=0.01
        converge=true;
        j=j+1;
    else
        mi=m;
        j=j+1;
    end
end
%}
%{
m1=0.87;
p1=0.4;
Py=955; %kips/K
Mp=4840; %Kin
%}
%{
syms m real

if p1>0
y1=(((1-m.^2)./(1+3.5*m.^2)).^0.5);
else
y1=-(((1-m.^2)./(1+3.5*m.^2)).^0.5); 
end

y2=(Mp*(m1-m).*(1+3.5*m.^2).*(2+7*m.^2));
y3=Py*(2*m.*(1+3.5*m.^2) + 7*m.*(1-m.^2));
yfinal=y1.*(1 + y2./y3) - p1;
[mcorrec]=vpasolve(yfinal==0,m,[0 1]);
if p1>0
   pcorrec=sqrt((1-mcorrec^2)/(1+3.5*mcorrec^2));
else
   pcorrec=-sqrt((1-mcorrec^2)/(1+3.5*mcorrec^2)); 
end
%mcorrec;
%pcorrec;
phi=mcorrec^2 + pcorrec^2 + 3.5*(mcorrec*pcorrec)^2;
%}
      
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
