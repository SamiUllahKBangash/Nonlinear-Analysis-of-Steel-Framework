function[lamda,buckledconfig]=lamdainelastic_rk2(lamdabar,n,Nel,Nnodes,elems,nodes,doffree,du,u,Pref,tol)
%Here RK-2 type incremental single step procedure is employed to evaluate
%secant Ke


dlamda=lamdabar/n;
fint=zeros(6,Nel);
fint_predictor=zeros(6,Nel);
for i=1:n    
    %-------Assembling Global Elastic Stiffness Matrix---------------------
    Ke=zeros(3*Nnodes);
    for iel=1:Nel
        elnodes=elems(iel,1:2);
        nodexy=nodes(elnodes,:); %nodal coordinates prior to load step
        eldofs=3*elnodes(1)-2:3*elnodes(1);
        eldofs=[eldofs 3*elnodes(2)-2:3*elnodes(2)];
        E=elems(iel,3);I=elems(iel,4); A=elems(iel,5);Fy=elems(iel,6);
        %stress_axial=(fint(1,iel)/A);
        if fint(1,iel) > 0
            stress_axial=(fint(1,iel)/A);
            if stress_axial>= 0.5*Fy && stress_axial<Fy
                E=4*E*((stress_axial/Fy)*(1 - (stress_axial/Fy)));
            elseif stress_axial>=Fy
                E=eps*E;
                %E=0;
            end
        end
        [ke_global,ke_local,T]=kelemelastic(nodexy,E,I,A);
        Ke(eldofs,eldofs)=Ke(eldofs,eldofs) + ke_global;
    end
    Keff=Ke(doffree,doffree);
    L=chol(Keff,'lower');       %Ke=LL' (Cholesky Factorization)
    Linv=L^-1; %eye(doffree);
    du(doffree)=(Linv'*Linv)*0.5*dlamda*(Pref(doffree));
    
    %----------Predictor Stage/LOCAL STATE DETERMINATION---------------------
    dfint=zeros(6,Nel);
    for iel=1:Nel
        elnodes=elems(iel,1:2);
        nodexy=nodes(elnodes,:); %System Current coordinates prior to load step.
        eldofs=3*elnodes(1)-2:3*elnodes(1);
        eldofs=[eldofs 3*elnodes(2)-2:3*elnodes(2)];
        E=elems(iel,3);I=elems(iel,4); A=elems(iel,5);
        if fint(1,iel) > 0
            stress_axial=(fint(1,iel)/A);
            if stress_axial>= 0.5*Fy && stress_axial<Fy
                E=4*E*((stress_axial/Fy)*(1 - (stress_axial/Fy)));
            elseif stress_axial>=Fy
                E=eps*E;
                %E=0;
            end
        end
        [ke_global,ke_local,T]=kelemelastic(nodexy,E,I,A);
        dfint(:,iel)=ke_local*T*du(eldofs);
        fint_predictor(:,iel)=fint(:,iel) + dfint(:,iel);
    end
    
    %----------Secant Stiffness Ke.Corrector Stage-----------------
    Ke_corrector=zeros(3*Nnodes);
    for iel=1:Nel
        elnodes=elems(iel,1:2);
        nodexy=nodes(elnodes,:); %nodal coordinates prior to load step
        eldofs=3*elnodes(1)-2:3*elnodes(1);
        eldofs=[eldofs 3*elnodes(2)-2:3*elnodes(2)];
        E=elems(iel,3);I=elems(iel,4); A=elems(iel,5);Fy=elems(iel,6);
        %stress_axial=(fint(1,iel)/A);
        if fint_predictor(1,iel) > 0
            stress_axial=(fint_predictor(1,iel)/A);
            if stress_axial>= 0.5*Fy && stress_axial<Fy
                E=4*E*((stress_axial/Fy)*(1 - (stress_axial/Fy)));
            elseif stress_axial>=Fy
                E=eps*E;
                %E=0;
            end
        end
        [ke_global,ke_local,T]=kelemelastic(nodexy,E,I,A);
        Ke_corrector(eldofs,eldofs)=Ke_corrector(eldofs,eldofs) + ke_global;
    end
    Keff=Ke_corrector(doffree,doffree);
    L=chol(Keff,'lower');       %Ke=LL' (Cholesky Factorization)
    Linv=L^-1; %eye(doffree);
    du(doffree)=(Linv'*Linv)*dlamda*(Pref(doffree));
    
    %----------Corrector Stage/LOCAL STATE DETERMINATION---------------------
    dfint=zeros(6,Nel);
    for iel=1:Nel
        elnodes=elems(iel,1:2);
        nodexy=nodes(elnodes,:); %System Current coordinates prior to load step.
        eldofs=3*elnodes(1)-2:3*elnodes(1);
        eldofs=[eldofs 3*elnodes(2)-2:3*elnodes(2)];
        E=elems(iel,3);I=elems(iel,4); A=elems(iel,5);
        if fint_predictor(1,iel) > 0
            stress_axial=(fint_predictor(1,iel)/A);
            if stress_axial>= 0.5*Fy && stress_axial<Fy
                E=4*E*((stress_axial/Fy)*(1 - (stress_axial/Fy)));
            elseif stress_axial>=Fy
                E=eps*E;
                %E=0;
            end
        end
        [ke_global,ke_local,T]=kelemelastic(nodexy,E,I,A);
        dfint(:,iel)=ke_local*T*du(eldofs);
        fint(:,iel)=fint(:,iel) + dfint(:,iel);
    end    
    
    u=u+du;
end

%------Kg Matrix Creation--------------
Kg=zeros(3*Nnodes);
for iel=1:Nel
   elnodes=elems(iel,1:2);
    nodexy=nodes(elnodes,:); %System Current coordinates prior to load step.    
    eldofs=3*elnodes(1)-2:3*elnodes(1);
    eldofs=[eldofs 3*elnodes(2)-2:3*elnodes(2)];
    Pj=fint(4,iel);
    [kg_global,kg_local,T]=kelem_geom(nodexy,Pj);
    Kg(eldofs,eldofs)=Kg(eldofs,eldofs) + kg_global;
end


Kgff=Kg(doffree,doffree);

%Conversion to Standard Form Hv=w*v
H=Linv*(-Kgff)*Linv'; % w=eigenvalue=1/lamda , v=eigenvector=L'*buckledconfig;

%Non-symmetric Form 
H2=(Linv'*Linv)*(-Kgff);

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
        buckledconfig=Linv'*y;
    else
        wo=w;
        yo=y;
        i=i+1;
    end
end
lamda=lamdacritical;
end
