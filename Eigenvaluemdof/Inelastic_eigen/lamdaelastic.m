function [lamda,fint]=lamdaelastic(Nel,Nnodes,elems,nodes,doffree,u,Pref,tol)
%{


%Welcome to Elastic Eigenvalue analysis for 2D Steel Frame Structures.
%Standard stiffness based finite element concepts are employed in a small-
%strain but moderate rotation setting. Cubic Hermitian Geometric Stiffness
%Matrix is used to approximate second order axial flexural coupling. Power
%Method, a type of vector iteration method, is applied to solve the
%standard Eigenvalue problem cast in Standard Form. (Ziemian 12.40).

%clc; clear all;

%---------------------------- Input File-----------------------------

[nodes,elems,bcs,loads]=example10point5;%getting input files
%elems array format:Node-near, Node-far, E, I, A, Fy, Zx

%---------------------System Variable Initialization-----------------------
for i=1:1
    Nel=size(elems,1);          %No of frame elements
    Nnodes=size(nodes,1);       %No of Nodes in the structure
    alldofs=1:3*Nnodes;
    
    %Nodes_updated=nodes;       %updated lagrangian formulation
    
    fint=zeros(6,Nel);          %trial/updated internal forces after load step
    %FINT=zeros(6,Nel);         %Converged Internal forces prior to load step
    %FINT_history=zeros(6,Nel);
    
    %du=zeros(3*Nnodes,1);      %incremental-displacement vector in global DSM coordinates
    u=zeros(3*Nnodes,1);        %displacement vector in global DSM coordinates
    %u_history=zeros(3*Nnodes,1);
    
    %zz=1;
    Pref=zeros(3*Nnodes,1);     %Proportional/Normalized force vector in Global DSM coordinates
    %lamda=0;
    %lamda_history(zz)=0;
end

%----------------Population of specified Pref DOFS--------------------------------
for ii=1:size(loads,1)
    thisdof=3*(loads(ii,1)-1)+loads(ii,2);
    Pref(thisdof)=loads(ii,3);
end

%---------------Population of contrained Displacement DOFS----------------

dofspec=[]; %constrained DOFs Array
for ii=1:size(bcs,1)
    thisdof=3*(bcs(ii,1)-1) + bcs(ii,2);
    dofspec=[dofspec thisdof];
    %u(thisdof)=bcs(ii,3);    %currently,the code doesn't account for nonlinear node settlement cases.
    %du(thisdof)=bcs(ii,3);
end
doffree=alldofs;
doffree(dofspec)=[]; %alldofs= (doffree + dofspec)


%}
%-------Assembling Global Elastic Stiffness Matrix---------------------
Ke=zeros(3*Nnodes);
for iel=1:Nel
    elnodes=elems(iel,1:2);
    nodexy=nodes(elnodes,:); %nodal coordinates prior to load step
    eldofs=3*elnodes(1)-2:3*elnodes(1);
    eldofs=[eldofs 3*elnodes(2)-2:3*elnodes(2)];
    E=elems(iel,3);I=elems(iel,4); A=elems(iel,5);    
    [ke_global,ke_local,T]=kelemelastic(nodexy,E,I,A);
    Ke(eldofs,eldofs)=Ke(eldofs,eldofs) + ke_global;
end

Keff=Ke(doffree,doffree);
L=chol(Keff,'lower');       %Ke=LL' (Cholesky Factorization)
Linv=L^-1; %eye(doffree);
u(doffree)=(Linv'*Linv)*(Pref(doffree));   

%----------LOCAL STATE DETERMINATION/Kg Matrix Creation---------------
fint=zeros(6,Nel);
Kg=zeros(3*Nnodes);
for iel=1:Nel
    elnodes=elems(iel,1:2);
    nodexy=nodes(elnodes,:); %System Current coordinates prior to load step.    
    eldofs=3*elnodes(1)-2:3*elnodes(1);
    eldofs=[eldofs 3*elnodes(2)-2:3*elnodes(2)];
    E=elems(iel,3);I=elems(iel,4); A=elems(iel,5);
    [ke_global,ke_local,T]=kelemelastic(nodexy,E,I,A);    
    fint(:,iel)=ke_local*T*u(eldofs);
    Pj=fint(4,iel);
    [kg_global,kg_local,T]=kelem_geom(nodexy,Pj);
    Kg(eldofs,eldofs)=Kg(eldofs,eldofs) + kg_global;
end

Kgff=Kg(doffree,doffree);

%Conversion to Standard Form Hv=w*v
H=Linv*(-Kgff)*Linv'; % w=eigenvalue=1/lamda , v=eigenvector=L'*buckledconfig;

%Non-symmetric Form 
%H2=(Linv'*Linv)*(-Kgff);

%Power Method Iterations
yo=ones(size(doffree,2),1);
%tol=1e-4;
converge=false;
wo=0;
i=1;
while ~converge
    yhat=H*yo;
    y=yhat/norm(yhat);
    w=y'*H*y;
    epsilon=abs((w-wo)/w)*100;
    if epsilon<=tol
        converge=true;
        lamdacritical=1/w;
        %buckledconfig=Linv'*y;
    else
        wo=w;
        yo=y;
        i=i+1;
    end
end

lamda=lamdacritical; %Critical Load Ratio Multiplier

end




    
    




