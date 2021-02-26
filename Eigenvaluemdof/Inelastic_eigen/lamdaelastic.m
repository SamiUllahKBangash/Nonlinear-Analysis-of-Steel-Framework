function [lamda,fint]=lamdaelastic(Nel,Nnodes,elems,nodes,doffree,u,Pref,tol)

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




    
    




