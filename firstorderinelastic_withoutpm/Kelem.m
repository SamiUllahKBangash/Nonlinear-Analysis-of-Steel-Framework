function [Kel] =  Kelem(nodexy, e, i, a,rho,mi,mj,mp)
%This function must return a 6x6 element stiffness matrix: [Kel]
%This matrix must be in Global Coordinates
%nodexy : [ x1  y1
%           x2  y2]
% E : Young's Modulus
% A : Area of cross-section
% I : Second moment of area

E1 = [ (nodexy(2,1)-nodexy(1,1))  (nodexy(2,2)-nodexy(1,2)) ];
L = norm(E1);
E1 = E1/L;
E2 = [-E1(2) E1(1)];
if abs(mp-abs(mi))<0.01*mp && abs(mp-abs(mj))>0.01*mp
%if abs(mi)>=0.99*mp && abs(mj)<0.99*mp    
    
     k_1=rho*e*i/L*[ a/i 0 0;     %perfectly elastic portion
                    0    4 2;
                    0    2 4];
     k_2=(1-rho)*e*i/L*[ a/i 0 0;  %EPP portion (subject to change)
                        0    0 0;
                        0    0 3];
    
elseif  abs(mp-abs(mi))>0.01*mp && abs(mp-abs(mj))<0.01*mp
%elseif abs(mi)<0.99*mp && abs(mj)>=0.99*mp    
     k_1=rho*e*i/L*[ a/i 0 0;     %perfectly elastic portion
                    0    4 2;
                    0    2 4];
     k_2=(1-rho)*e*i/L*[ a/i 0 0;  %EPP portion (subject to change)
                        0    3 0;
                        0    0 0];
elseif  abs(mp-abs(mi))<0.01*mp && abs(mp-abs(mj))<0.01*mp
%elseif abs(mi)>=0.99*mp && abs(mj)>=0.99*mp    
     k_1=rho*e*i/L*[ a/i 0 0;     %perfectly elastic portion
                    0    4 2;
                    0    2 4];
     k_2=(1-rho)*e*i/L*[ a/i 0 0;  %EPP portion (subject to change)
                        0    0 0;
                        0    0 0];
else
    
    k_1=rho*e*i/L*[ a/i 0 0;     %perfectly elastic portion
                    0    4 2;
                    0    2 4];
     k_2=(1-rho)*e*i/L*[ a/i 0 0;  %EPP portion (subject to change)
                        0    4 2;
                        0    2 4];
end
    
 %{    
Kel_LOC1 = [ A*E/L       0         0   -A*E/L         0              0;
            0       3*E*I/(L^3)   0     0       -3*E*I/(L^3)   3*E*I/(L^2); 
            0           0         0     0            0              0;
           -A*E/L       0         0   A*E/L          0              0;
            0      -3*E*I/(L^3)   0     0       3*E*I/(L^3)    -3*E*I/(L^2);
        	0       3*E*I/(L^2)   0     0      -3*E*I/(L^2)     3*E*I/L ];
        
Kel_LOC2 = 0.00001*[ A*E/L    0   0   -A*E/L  0   0;
        0   12*E*I/(L^3)    6*E*I/(L^2)   0    -12*E*I/(L^3)   6*E*I/(L^2); 
        0   6*E*I/(L^2)     4*E*I/L      0    -6*E*I/(L^2)    2*E*I/L;
        -A*E/L    0   0   A*E/L  0   0;
        0   -12*E*I/(L^3)   -6*E*I/(L^2)  0    12*E*I/(L^3)    -6*E*I/(L^2);
        0   6*E*I/(L^2)     2*E*I/L    0     -6*E*I/(L^2)    4*E*I/L ];
%}
k=k_1 + k_2;
at=[1  0   0 -1  0    0;             %element transformation matrix
    0 1/L 1 0 -1/L 0;
    0 1/L 0 0 -1/L 1];
Qrot = [E1 ; E2];
Qrot(3,3) = 1;
Tmatrix = [Qrot zeros(3); zeros(3) Qrot];
Kel = Tmatrix'*(at'*k*at)*Tmatrix;
end