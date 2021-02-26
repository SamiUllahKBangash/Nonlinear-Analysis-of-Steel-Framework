%Welcome to In-Elastic Eigenvalue analysis for 2D Steel Frame Structures.
%Standard stiffness based finite element concepts are employed in a small-
%strain but moderate rotation setting. Cubic Hermitian Geometric Stiffness
%Matrix is used to approximate second order axial flexural coupling. Power
%Method, a type of vector iteration method, is applied to solve the
%standard Eigenvalue problem cast in Standard Form. (Ziemian 12.40).

clc; clear all;

%---------------------------- Input File-----------------------------

[nodes,elems,bcs,loads]=example10point5;%getting input files
%elems array format:Node-near, Node-far, E, I, A, Fy, Zx

%--------Input Parameters for critical Load Ratio Optimization-----------
n=150         %no of incremental steps for Stiffness Matrix and fint Calculation
tol=1e-6       %tolerance parameter for Power Iteration Convergence
tau_factor=0.1 %lower bound lamdabar factor for falsi iteration
rk2_flag=true; %Flag to activate RK-2 type incremental Ke/Kg solver

%---------------------System Variable Initialization-----------------------
for i=1:1
    Nel=size(elems,1);          %No of frame elements
    Nnodes=size(nodes,1);       %No of Nodes in the structure
    alldofs=1:3*Nnodes;
    
    %Nodes_updated=nodes;       %updated lagrangian formulation
    
    %fint=zeros(6,Nel);         %trial/updated internal forces after load step
    FINT=zeros(6,Nel);          %Converged Internal forces prior to load step
    FINT_history=zeros(6,Nel);
    
    du=zeros(3*Nnodes,1);       %incremental-displacement vector in global DSM coordinates
    u=zeros(3*Nnodes,1);        %displacement vector in global DSM coordinates
    u_history=zeros(3*Nnodes,1);
    
    zz=1;
    Pref=zeros(3*Nnodes,1);     %Proportional/Normalized force vector in Global DSM coordinates
    %lamda=0;
    %lamda_history(zz)=0;
    
end

%----------------Population of Specified Pref DOFS--------------------------------
for ii=1:size(loads,1)
    thisdof=3*(loads(ii,1)-1)+loads(ii,2);
    Pref(thisdof)=loads(ii,3);
end

%---------------Population of Constrained Displacement DOFS----------------

dofspec=[]; %constrained DOFs Array
for ii=1:size(bcs,1)
    thisdof=3*(bcs(ii,1)-1) + bcs(ii,2);
    dofspec=[dofspec thisdof];
end
doffree=alldofs;
doffree(dofspec)=[]; %alldofs= (doffree + dofspec)


%-----Running Elastic Eigenvalue first--------
[lamda,fint]=lamdaelastic(Nel,Nnodes,elems,nodes,doffree,u,Pref,tol);

inelastic=[];
%i=1;
for iel=1:Nel
    Py=elems(iel,5)*elems(iel,6);
    if lamda*fint(1,iel)>0 && lamda*fint(1,iel)>= 0.5*Py
        inelastic=[inelastic iel];
    end
end

if isempty(inelastic)
    criticalratio_elastic=lamda
    disp('Elastic Buckling governs. Program will terminate')
    return;
else
    disp('Inelastic Buckling governs.')
end


criticalratio_elastic=lamda
lamdaElastic=lamda;

%-------Initiation of Regular Falsi iterations----------------
dlamda=tau_factor*lamdaElastic - lamdaElastic;
tau_a=0;
tau_b=1;
converge=false;
j=1;
while ~converge
    
    if ~rk2_flag
        lamdabar=lamdaElastic+(tau_a*dlamda);
        [lamda,~]=lamdainelastic(lamdabar,n,Nel,Nnodes,elems,nodes,doffree,du,u,Pref,tol);
        lamda_a=lamda;
        
        lamdabar=lamdaElastic+(tau_b*dlamda);
        [lamda,~]=lamdainelastic(lamdabar,n,Nel,Nnodes,elems,nodes,doffree,du,u,Pref,tol);
        lamda_b=lamda;
        
        tau_r= tau_b - ((lamda_b-1)*(tau_a-tau_b))/(lamda_a - lamda_b);
        
        lamdabar=lamdaElastic + (tau_r*dlamda);
        [lamda,buckledconfig]=lamdainelastic(lamdabar,n,Nel,Nnodes,elems,nodes,doffree,du,u,Pref,tol);
        
        lamda_r=lamda;
    
    else 
        lamdabar=lamdaElastic+(tau_a*dlamda);
        [lamda,~]=lamdainelastic_rk2(lamdabar,n,Nel,Nnodes,elems,nodes,doffree,du,u,Pref,tol);
        lamda_a=lamda;
        
        lamdabar=lamdaElastic+(tau_b*dlamda);
        [lamda,~]=lamdainelastic_rk2(lamdabar,n,Nel,Nnodes,elems,nodes,doffree,du,u,Pref,tol);
        lamda_b=lamda;
        
        tau_r= tau_b - ((lamda_b-1)*(tau_a-tau_b))/(lamda_a - lamda_b);
        
        lamdabar=lamdaElastic + (tau_r*dlamda);
        [lamda,buckledconfig]=lamdainelastic_rk2(lamdabar,n,Nel,Nnodes,elems,nodes,doffree,du,u,Pref,tol);
        
        lamda_r=lamda;
    end
        
    if abs(lamda_r-1)<=tol
        converge=true;
        break;
    end
    if (lamda_r-1)*(lamda_b-1)>0
        tau_b=tau_r;
    elseif (lamda_r-1)*(lamda_a-1)>0
        tau_a=tau_r;
    end
    j=j+1;
end

criticalratio_inelastic=lamdabar
buckledconfig























