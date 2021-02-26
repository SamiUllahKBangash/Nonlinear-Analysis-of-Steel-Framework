function [kg_global,kg_local,T]=kelem_geom(nodexy,Pj)
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

Qrot = [E1(1) E1(2) 0;
      -E1(2)  E1(1) 0;
        0       0   1];
T = [Qrot zeros(3); 
    zeros(3) Qrot];

        
kg=(Pj/L)*[1  0        0       -1   0      0;
            0  6/5     L/10      0  -6/5    L/10;
            0  L/10 (2*L^2)/15   0  -L/10 -L^2/30;
           -1   0       0        1    0      0;
            0  -6/5   -L/10      0  6/5   -L/10;
            0  L/10  -L^2/30     0 -L/10  (2*L^2)/15];
 
        
kg_local=kg;
kg_global=T'*kg_local*T;
end


    
        

