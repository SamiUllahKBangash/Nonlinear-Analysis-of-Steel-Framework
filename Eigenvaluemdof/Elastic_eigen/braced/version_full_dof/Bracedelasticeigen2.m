%Welcome. This routine performs elastic critical buckling load analysis of
%UNBRACED-Planar-Steel-Frames. Multiple options are provided to run critical
%analysis namely Power Method, polynomial expansion and Stability
%Functions. Power Method and polynomial expansion require the use of Cubic
%Hermitian Matrix to approximate Geometric Stiffness effect. Overall goal
%is to compare the performance of Cubic Hermitian Matric with the more
%accurate Stability Functions. Also, the most important benefit of Power
%Method is highlighted i.e. its ability to handle arbitrarily large frame
%sizes due to superior performance/stability and convergence rate

%Note: This is the full version of Braced Frame Elastic Critical Analysis 
%Routine. Here the standard DOF numbering of 3DOF per node is retained.
%Consequently, the standard geometric and stability matrices are employed
%and bracing condition is defined explicitly in the input file. 


%clc; clear all;

%---------------------------- Input File-----------------------------
[nodes,StoreyNodes,elems,bcs,~,loads]=example10point5;
%[nodes,elems,bcs,loads]=Fullonebayframe;%getting input files
%elems array format:Node-near, Node-far, E, I, A, Fy, Zx
flag_stability=true;
if ~flag_stability
    flag_polyexpansion=true;
end
%---------------------System Variable Initialization-----------------------
for i=1:1
    Nel=size(elems,1);          %No of frame elements
    Nnodes=size(nodes,1);       %No of Nodes in the structure
    alldofs=1:3*Nnodes;
    fint=zeros(6,Nel);          %trial/updated internal forces after load step
    u=zeros(3*Nnodes,1);        %displacement vector in global DSM coordinates
    Pref=zeros(3*Nnodes,1);     %Proportional/Normalized force vector in Global DSM coordinates
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
end
doffree=alldofs;
doffree(dofspec)=[]; %alldofs= (doffree + dofspec)

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
L=chol(Keff,'lower');
Linv=L^-1; %eye(doffree);
u(doffree)=Linv'*Linv*(Pref(doffree));


if ~flag_stability
    Kg=zeros(3*Nnodes);
else
    Kg=sym(zeros(3*Nnodes));
    syms x real
end

%-----Local State Determination, Kg,Kt formation--------------

for iel=1:Nel
    elnodes=elems(iel,1:2);
    nodexy=nodes(elnodes,:);
    eldofs=3*elnodes(1)-2:3*elnodes(1);
    eldofs=[eldofs 3*elnodes(2)-2:3*elnodes(2)];
    E=elems(iel,3);I=elems(iel,4); A=elems(iel,5);
    [ke_global,ke_local,T]=kelemelastic(nodexy,E,I,A);
    fint(:,iel)=ke_local*T*u(eldofs);
    Pj=fint(4,iel);
    if ~flag_stability
        [kg_global,kg_local,~]=kelem_geom(nodexy,Pj);
    else
        [kg_global,kg_local,T]=kelem_geom_exact(nodexy,Pj,E,A,I);
        
    end
    
    Kg(eldofs,eldofs)=Kg(eldofs,eldofs) + kg_global;
end

Kgff=Kg(doffree,doffree);

if ~flag_stability
    %Conversion to Standard Form Hv=w*v
    H=Linv*(-Kgff)*Linv'; % w=eigenvalue=1/lamda , v=eigenvector=Linv'*buckledconfig;
    
    %Non-symmetric Form
    %H2=(Linv'*Linv)*(-Kgff);
    
    %Power Method Iterations
    
    if ~flag_polyexpansion
        %power Method Iterations
        yo=ones(size(doffree,2),1);
        tol=1e-4;
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
                buckledconfig=Linv'*y;
            else
                wo=w;
                yo=y;
                i=i+1;
            end
        end
    else %Polynomial expansion of standard form |H-wI|=0
        syms w real
        z=size(H,1);
        eq1=det(H-w*eye(z));
        sol=vpasolve(eq1);
        lamdacritical=1/(max(sol(sol>0)));
    end
    
else
    eq1=det(Kgff);
    %eq1=@(y)subs(eq1,x,y);
    sol=[];
    for i=1:4
        sol=[sol vpasolve(eq1,x,'random',true)];
        %sol=vpasolve(eq1,x);
        %sol=fzero(eq1,[1930 7702]);
    end
    lamdacritical=min(sol(sol>0));
    %lamdacritical=sol;
end
lamdacritical
%buckledconfig











