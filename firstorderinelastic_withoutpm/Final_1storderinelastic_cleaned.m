clear all; clc; close all; %clear all existing variables (new start)

%Get the input file
filename=uigetfile('*.txt','Select input file');
[nodes, elems, bcs, loads, E, I, A, Mp] = FOIgetframedata2D(filename);  %rho is fixed (0.00001)


Nel = size(elems,1);
Nnodes = size(nodes,1);

%Decide degrees of freedom + Initialize Matrices
alldofs=1:3*Nnodes;


%Note: DOFs corresponding to node "i" are
% [3*(i-1)+1   3*(i-1)+2    3*(i-1)+3]
duini = zeros(3*Nnodes,1);
uini = zeros(3*Nnodes,1);
fini = zeros(3*Nnodes,1);

%Nodal Loads
dofload=[];
for ii = 1:size(loads,1)
    thisdof=3*(loads(ii,1)-1)+loads(ii,2);
    fini(thisdof) = loads(ii,3);
    dofload=[dofload thisdof];
end

%Boundary Conditions

dofspec = [];
for ii = 1:size(bcs,1)
    thisdof = 3*(bcs(ii,1)-1)+bcs(ii,2);
    dofspec = [dofspec thisdof];
    uini(thisdof) = bcs(ii,3);
end
doffree = alldofs;
doffree(dofspec) = []; %Delete sepcified dofs from All dofs

rho=0.0001; %hardening parameter 
zz=1; %counter to update situation
Mint=zeros(2,Nel); %column shows element no and row shows near and far end end moment values
disp_node2x=0; %user defines DOF no for tracking its 1st order inelastic path
H=0;
H_history(zz)=0;
lateral_disp_node2x(zz)=0;
Mint_history=[];
SF_history=[];
SF_govern_history=[];

Collapse=0; %flag to check stiffness condition number
while Collapse==0
%for iter=1:3
K = zeros(3*Nnodes);
%f=zeros(3*Nnodes,1);

%Initialize Global Stiffness Matrix
for iel = 1:Nel
    elnodes = elems( iel, 1:2);
    nodexy = nodes(elnodes, :);
    e=E(iel); i=I(iel); a=A(iel);
    mi=Mint(1,iel); mj=Mint(2,iel); mp=Mp(iel);
    [Kel] = Kelem(nodexy, e, i, a,rho,mi,mj,mp);
          
    %Assemble element stiffness matrix into Global stiffness matrix K
    eldofs = 3*(elnodes(1)-1)+1:3*(elnodes(1));
    eldofs = [eldofs 3*(elnodes(2)-1)+1:3*elnodes(2)];
    K(eldofs,eldofs) = K(eldofs,eldofs) + Kel;
end

Keff=K(doffree,doffree);

%Solve
Dispff = Keff\fini(doffree);      %Incremental Displacements for free dofs
duini(doffree) = Dispff;
dMint=zeros(2,Nel);
sf=zeros(2,Nel);



for iel = 1:Nel
    elnodes = elems( iel, 1:2);             %To find elemental forces based on these displacements
    nodexy = nodes(elnodes, :);
    e=E(iel); i=I(iel); a=A(iel);
     mi=Mint(1,iel); mj=Mint(2,iel); mp=Mp(iel);
    [Kel] = Kelem(nodexy, e, i, a,rho,mi,mj,mp);
    

    ReqDef = [duini((elnodes(1)*3)-2); duini((elnodes(1)*3)-1); duini(elnodes(1)*3); duini((elnodes(2)*3)-2); duini((elnodes(2)*3)-1); duini(elnodes(2)*3)];
    dMint(1,iel)= Kel(3,:)*ReqDef;
    dMint(2,iel)=Kel(6,:)*ReqDef;
     
    for j=1:2 
                     %if Mint(j,i)/elem_moment(j,i) < 0
                     %print('elastic unloading occuring!')
                     %  return
           if abs(dMint(j,iel))>0.01
              sf(j,iel)=(Mp(iel)-abs(Mint(j,iel)))/abs(dMint(j,iel));
           end
    end
    
end
SF_history=[SF_history;sf];
zz=zz+1;
 %scale factor that will be used to scale up the current reference load, nodal displacements,elem_moment
 SF_govern=min(sf(sf > 0));
 SF_govern_history=[SF_govern_history SF_govern];
 Mint=Mint+SF_govern*(dMint);
 Mint_history=[Mint_history;Mint];
 disp_node2x=disp_node2x+SF_govern*duini(4);
 H=H+SF_govern*fini(4);
 H_history(zz)=H;
 lateral_disp_node2x(zz)=disp_node2x;
 
 if (lateral_disp_node2x(zz)/lateral_disp_node2x(zz-1)) >=100 && zz>2
     disp 'collapse imminent'
     Collapse=1;
 end

end
 


