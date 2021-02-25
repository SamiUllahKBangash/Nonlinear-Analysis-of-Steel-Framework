function [kt,T,G,ke,km]=kelemp(nodexy,E,I,A,Py,Mp,Pi,Mi,Pj,Mj,tol)
%This function must return a 6x6 element tangent stiffness matrix: [kt]
%This matrix must be in Global Coordinates
%nodexy : [ x1  y1
%           x2  y2]
% E : Young's Modulus
% A : Area of cross-section
% I : Second moment of area

E1 =nodexy(2,:) - nodexy(1,:);
L = norm(E1);
E1 = E1/L;
E2 = [-E1(2) E1(1)];
Qrot = [E1 ; E2];
Qrot(3,3) = 1;
T = [Qrot zeros(3); zeros(3) Qrot];
pi=Pi/Py; mi=Mi/Mp; pj=Pj/Py; mj=Mj/Mp;
phi_i=pi^2 + mi^2 + 3.5*(pi*mi)^2;
phi_j=pj^2 + mj^2 + 3.5*(pj*mj)^2;
dphipi=(2*pi/Py)+(7*pi*mi^2/Py);
dphimi=(2*mi/Mp)+(7*pi^2*mi/Mp);
dphipj=(2*pj/Py)+(7*pj*mj^2/Py);
dphimj=(2*mj/Mp)+(7*pj^2*mj/Mp);

ke =      [ A*E/L            0               0       -A*E/L       0             0;
            0          12*E*I/(L^3)    6*E*I/(L^2)   0    -12*E*I/(L^3)   6*E*I/(L^2); 
            0           6*E*I/(L^2)     4*E*I/L      0     -6*E*I/(L^2)    2*E*I/L;
          -A*E/L            0             0         A*E/L       0              0;
            0           -12*E*I/(L^3)   -6*E*I/(L^2)  0    12*E*I/(L^3)    -6*E*I/(L^2);
            0            6*E*I/(L^2)     2*E*I/L    0     -6*E*I/(L^2)      4*E*I/L ];

if phi_i<(1-tol) && phi_j<(1-tol) 
    G=0;
    km=zeros(6);
    
elseif abs(phi_i-1)<=tol && phi_j<(1-tol)    
    G=[dphipi;0;dphimi;0;0;0];
    km=-ke*G*((G'*ke*G)^-1)*G'*ke;
    
elseif phi_i<(1-tol) && abs(phi_j-1)<=tol 
    G=[0;0;0;dphipj;0;dphimj;];
    km=-ke*G*((G'*ke*G)^-1)*G'*ke;
    
elseif abs(phi_i-1)<=tol && abs(phi_j-1)<=tol 
    G=[dphipi 0;
        0     0;
       dphimi 0;
        0 dphipj;
        0     0;
        0 dphimj];         
    km=-ke*G*((G'*ke*G)^-1)*G'*ke;
end

kt= T'*(ke+km)*T;

end