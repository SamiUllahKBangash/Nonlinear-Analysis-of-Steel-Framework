function [ke_global,ke_local,T]=kelemelastic(nodexy,E,I,A)
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
%E2 = [-E1(2) E1(1)];
Qrot = [E1(1) E1(2) 0;
      -E1(2)  E1(1) 0;
        0       0   1];
%Qrot(3,3) = 1;
T = [Qrot zeros(3); 
    zeros(3) Qrot];

ke =    [ A*E/L            0               0       -A*E/L       0             0;
            0          12*E*I/(L^3)    6*E*I/(L^2)   0    -12*E*I/(L^3)   6*E*I/(L^2); 
            0           6*E*I/(L^2)     4*E*I/L      0     -6*E*I/(L^2)    2*E*I/L;
          -A*E/L            0             0         A*E/L       0              0;
            0           -12*E*I/(L^3)   -6*E*I/(L^2)  0    12*E*I/(L^3)    -6*E*I/(L^2);
            0            6*E*I/(L^2)     2*E*I/L    0     -6*E*I/(L^2)      4*E*I/L ];
        

ke_local=ke;
ke_global=T'*ke_local*T;
end


    
        

