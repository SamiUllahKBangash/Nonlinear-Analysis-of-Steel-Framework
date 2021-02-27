function [ke_global,ke_local,T]=Braced_kelemelastic(nodexy,E,I)
%This function must return a 6x6 element tangent stiffness matrix: [kt]
%This matrix must be in Global Coordinates
%nodexy : [ x1  y1
%           x2  y2]
% E : Young's Modulus
% A : Area of cross-section
% I : Second moment of area
E1 =nodexy(2,:) - nodexy(1,:);
L = norm(E1);
T=eye(2);
ke =  E*I/L*[4  2;            
             2  4];
        
ke_local=ke;
ke_global=T'*ke_local*T;
end
