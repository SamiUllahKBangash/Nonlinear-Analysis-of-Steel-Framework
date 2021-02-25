function [kg_global,kg_local,T]=kelem_geom_exact(nodexy,Pj,E,I)
%This function must return a 6x6 element tangent stiffness matrix: [kt]
%This matrix must be in Global Coordinates
%nodexy : [ x1  y1
%           x2  y2]
% E : Young's Modulus
% A : Area of cross-section
% I : Second moment of area
E1 =nodexy(2,:) - nodexy(1,:);
L = norm(E1);
%E1 = E1/L;
%E2 = [-E1(2) E1(1)];
%Qrot = [E1(1) E1(2) 0;
 %   -E1(2)  E1(1) 0;
  %  0       0   1];
%Qrot(3,3) = 1;
%T = [Qrot zeros(3);
 %   zeros(3) Qrot];
T=eye(2);

syms x real
tol=1e-6;
%tol=0;
Pe=pi^2*E*I/L^2;
P=x*abs(Pj);
lamda=pi*sqrt(P/Pe);

if Pj<(0-tol)
    phi_c=2-2*cos(lamda)-lamda*sin(lamda);
    %phi_1=lamda^3*sin(lamda)/(12*phi_c);
    %phi_2=lamda^2*(1-cos(lamda))/(6*phi_c);
    phi_3=lamda*(sin(lamda) - lamda*cos(lamda))/(4*phi_c);
    phi_4=lamda*(lamda - sin(lamda))/(2*phi_c);
    
elseif Pj>(0+tol)
    phi_t=2 - 2*cosh(lamda) + lamda*sinh(lamda);
    %phi_1=lamda^3*sinh(lamda)/(12*phi_t);
    %phi_2=lamda^2*(cosh(lamda)-1)/(6*phi_t);
    phi_3=lamda*(lamda*cosh(lamda) -sinh(lamda))/(4*phi_t);
    phi_4=lamda*(sinh(lamda) - lamda)/(2*phi_t);
else
    %phi_1=1;%phi_2=1;
    phi_3=1;phi_4=1;
end


kg=(E*I/L)*[4*phi_3   2*phi_4;   
            2*phi_4   4*phi_3];

kg_local=kg;
kg_global=T'*kg_local*T;
end






