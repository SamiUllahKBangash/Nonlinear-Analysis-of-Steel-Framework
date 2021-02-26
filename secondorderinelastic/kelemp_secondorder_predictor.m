function[kt_G,kt_L,ke_L,kg_L,km_L,T,G]=kelemp_secondorder_predictor(nodexy,E,I,A,Py,Mp,Forces_internal_prev,Forces_internal_midstep,tol)

%This function must return a 6x6 element SECANT (RK-2) stiffness matrix:[ksec]
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


ke_L =    [ A*E/L            0               0       -A*E/L       0             0;
    0          12*E*I/(L^3)    6*E*I/(L^2)   0    -12*E*I/(L^3)   6*E*I/(L^2);
    0           6*E*I/(L^2)     4*E*I/L      0     -6*E*I/(L^2)    2*E*I/L;
    -A*E/L            0             0         A*E/L       0              0;
    0           -12*E*I/(L^3)   -6*E*I/(L^2)  0    12*E*I/(L^3)    -6*E*I/(L^2);
    0            6*E*I/(L^2)     2*E*I/L    0     -6*E*I/(L^2)      4*E*I/L ];
kg_L=(Forces_internal_midstep(3)/L)*[1  0        0       -1   0      0;
    0  6/5     L/10      0  -6/5    L/10;
    0  L/10 (2*L^2)/15   0  -L/10 -L^2/30;
    -1   0       0        1    0      0;
    0  -6/5   -L/10      0  6/5   -L/10;
    0  L/10  -L^2/30     0 -L/10  (2*L^2)/15];

pi_prev=Forces_internal_prev(1)/Py;
mi_prev=Forces_internal_prev(2)/Mp;
pj_prev=Forces_internal_prev(3)/Py;
mj_prev=Forces_internal_prev(4)/Mp;

pi_mid=Forces_internal_midstep(1)/Py;
mi_mid=Forces_internal_midstep(2)/Mp;
pj_mid=Forces_internal_midstep(3)/Py;
mj_mid=Forces_internal_midstep(4)/Mp;

phi_i_prev=pi_prev^2 + mi_prev^2 + 3.5*(pi_prev*mi_prev)^2;
phi_j_prev=pj_prev^2 + mj_prev^2 + 3.5*(pj_prev*mj_prev)^2;

if abs(phi_i_prev -1)<=tol
    pval=pi_mid;
    mval=mi_mid;
    [tau_r]=falsi_midstep_surface(pval,mval,tol);
    dphipi=(1/(1+tau_r)^2)*((2*pi_mid/Py)   +   (7/(1+tau_r)^2)*(pi_mid*mi_mid^2/Py));
    dphimi=(1/(1+tau_r)^2)*((2*mi_mid/Mp)   +   (7/(1+tau_r)^2)*(pi_mid^2*mi_mid/Mp));
end

if abs(phi_j_prev -1)<=tol
    pval=pj_mid;
    mval=mj_mid;
    [tau_r]=falsi_midstep_surface(pval,mval,tol);
    dphipj=(1/(1+tau_r)^2)*((2*pj_mid/Py)   +   (1/(1+tau_r)^2)*(7*pj_mid*mj_mid^2/Py));
    dphimj=(1/(1+tau_r)^2)*((2*mj_mid/Mp)   +   (1/(1+tau_r)^2)*(7*pj_mid^2*mj_mid/Mp));
end




if phi_i_prev<(1-tol) && phi_j_prev<(1-tol)
    G=0;
    km_L=zeros(6);
elseif abs(phi_i_prev-1)<=tol && phi_j_prev<(1-tol)
    G=[dphipi;0;dphimi;0;0;0];
    km_L=-(ke_L+kg_L)*G*((G'*(ke_L+kg_L)*G)^-1)*G'*(ke_L+kg_L);
elseif phi_i_prev<(1-tol) && abs(phi_j_prev-1)<=tol
    G=[0;0;0;dphipj;0;dphimj;];
    km_L=-(ke_L+kg_L)*G*((G'*(ke_L+kg_L)*G)^-1)*G'*(ke_L+kg_L);
elseif abs(phi_i_prev-1)<=tol && abs(phi_j_prev-1)<=tol
    G=[dphipi 0;
        0     0;
        dphimi 0;
        0 dphipj;
        0     0;
        0 dphimj];
    km_L=-(ke_L+kg_L)*G*((G'*(ke_L+kg_L)*G)^-1)*G'*(ke_L+kg_L);
end
kt_L=ke_L+kg_L+km_L;
kt_G= T'*kt_L*T;
end
