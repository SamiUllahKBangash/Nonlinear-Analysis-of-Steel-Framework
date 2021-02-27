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

T=eye(2);

 kg=(Pj/L)*[(2*L^2)/15   -L^2/30;           
            -L^2/30   (2*L^2)/15];   
    
        
kg_local=kg;
kg_global=T'*kg_local*T;
end
