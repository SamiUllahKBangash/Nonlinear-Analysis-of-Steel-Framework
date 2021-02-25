%eq1=@(x)(x*(sin(x) - x*cos(x)))/(2- 2*cos(x) - x*sin(x));
%eq2=@(x)(2 - ((x-sin(x))/(sin(x)-x*cos(x)))^2);
%eq3=@(x)eq1*eq2;
%eq2=@(x)((x*(sin(x) - x*cos(x)))/(2- 2*cos(x) - x*sin(x)))^2;
%x0=[0.01 2*pi];


%{
%----code for the stick model------
syms x real
s=x*(sin(x) - x*cos(x))/(2 - 2*cos(x) - x*sin(x));
c=(x - sin(x))/(sin(x) - x*cos(x));
K=[2*s s*c;
   s*c s];
%det(K)
eq1=det(K);

sol=[];
for i=1:10
sol=[sol vpasolve(eq1,x,[0 2*pi],'random',true)];
end
solmin=min(sol)
%eq1=@(x)(2*x^3*sin(x) + 2*x^4*cos(x)^2 + x^2*sin(x)^2 - x^4 - 4*x^3*cos(x)*sin(x))/(2*cos(x) + x*sin(x) - 2)^2;

%eq1=(x*(sin(x) - x*cos(x)))/(2- 2*cos(x) - x*sin(x))*(2 - ((x-sin(x))/(sin(x)-x*cos(x)))^2);
%eq1=@(x)(x*(sin(x) - x*cos(x)))/(2- 2*cos(x) - x*sin(x))*(2 - ((x-sin(x))/(sin(x)-x*cos(x)))^2);
x0=pi;
fzero(eq1,x0)
%}

%------Sway Frame model----------
%{
syms x real
sii=(x*sin(x) - x^2*cos(x))/(2 - 2*cos(x) - x*sin(x));
sij=(x^2 - x*sin(x))/(2 - 2*cos(x) - x*sin(x));
S=sii+sij;
K=[-(sii+6) S;
    S x^2-2*S];
eq1=det(K);
sol=[];
for i=1:10
sol=[sol vpasolve(eq1,x,[pi/2 pi],'random',true)];
end
solmin=min(sol)
%}
syms x real
sii=(x*sin(x) - x^2*cos(x))/(2 - 2*cos(x) - x*sin(x));
sij=(x^2 - x*sin(x))/(2 - 2*cos(x) - x*sin(x));
eq1=sii - (sij^2/sii) + 2;
sol=[];
for i=1:2
sol=[sol vpasolve(eq1,x,[pi pi/0.7],'random',true)];
end
solmin=min(sol(sol>0))





