%Welcome. In this third take on first order inelastic analysis, we will be
%updating system variables without re-running the entire Tangent stiffness
%matrix assembly process all over again. i.e. after computing critical tau factor,
%the existing dfint values are scaled by tau directly instead of repeating
%the same load step (with tau correction) from while loop. After scaling
%dfint is done, we update fint_trial,phi arrays, lamdas,and displacement
%arrays by tau as well and move to next load step with some user defined
%dlamda_1 (subject to user discretion.)




clc;clear all;
%format long
%--------------------------Input file-------------------------------------
[nodes,elems,bcs,loads]=example10point5;
%elems array format:Node-near, Node-far, E, I, A, Fy, Zx

%------------------User defined Analysis Parameters-----------------------
disp('User defined Analysis Parameters')
lamdamax=350                       %maximum applied load ratio assigned by the analyst.
dlamda_1=50                        %incremental load ratio/step. |dlamda_1| must be < than lamdamax
tol=1e-4                           %tolerance on yield surface convergence
tol2=1e-8                          %tolerance for point projection on yield surface: only applicable to radial point projection algorithm
point_projection_radial=true       %radial point projection w.r.t. pm curve origin
point_projection_normal=false      %orthogonal point projection along surface gradient.
%tol2 defined above must be <= tol


%--------------------System Variable Initialization------------------------
for i=1:1
    Nel=size(elems,1);             %No of frame elements
    Nnodes=size(nodes,1);          %No of Nodes in the structure
    alldofs=1:3*Nnodes;
    
    fint=zeros(6,Nel);             %trial internal forces after load step
    FINT=zeros(6,Nel);             %storing converged internal forces prior to loadstep
    FINT_history=zeros(6,Nel);     %storing history of converged internal forces
    
    du=zeros(3*Nnodes,1);
    u=zeros(3*Nnodes,1);           %total nodal displacements/rotations in DSM coordinates
    u_history=zeros(3*Nnodes,1);
    
    up=zeros(6,Nel);               %plastic member-end displacements/rotations (in local axis)
    up_history=zeros(6,Nel);
    u_total=zeros(6,Nel);          %total member-end displacements/rotations (in local axis)
    u_total_history=zeros(6,Nel);
    
    n=0;
    zz=1;
    Pref=zeros(3*Nnodes,1);
    lamda=0;
    lamda_history(zz)=0;
    
    phi=zeros(6,Nel);              %array to store yielding condition of both member ends.
    %phi:[currentphi-near,previousphi-near,currentphi_tol-near,currentphi-far,previousphi-far,currentphi_tol-far)]'
    phi_history=zeros(6,Nel);
    
    hingetable=[];                 %bookkeeping of hinges as they form. Format: element# node # lamdaval
end
%----------------Population of Specified Pref DOFS-------------------------
dofload=[];
for ii=1:size(loads,1)
    thisdof=3*(loads(ii,1)-1)+loads(ii,2);
    Pref(thisdof)=loads(ii,3);   %normalized proportional loads
end

%------------------Population of Constrained Displacement DOFS------------
dofspec=[];
for ii=1:size(bcs,1)
    thisdof=3*(bcs(ii,1)-1)+bcs(ii,2);
    dofspec=[dofspec thisdof];
    %u(thisdof)=bcs(ii,3);     %currently,the code doesn't account for nonlinear node settlement cases.
    %du(thisdof)=bcs(ii,3);
end
doffree=alldofs;
doffree(dofspec)=[];           %alldofs = doffree + dofspec



%------------------Main Program Load Increment Loop-----------------------
while lamda<=lamdamax
    
    %initializing load-step magnitude
    if lamda==lamdamax
        disp('Maximum Applied Load Ratio reached. Program will terminate')
        break;
    end
    if (lamda+dlamda_1)>lamdamax
        dlamda=lamdamax-lamda;
    else
        dlamda=dlamda_1;
    end
    
    %------Creation of Tangent Stiffness Matrix at the start of load step--
    Kt=zeros(3*Nnodes);
    for iel=1:Nel
        elnodes=elems(iel,1:2);
        nodexy=nodes(elnodes,:);
        eldofs=3*(elnodes(1)-1)+1:3*elnodes(1);
        eldofs=[eldofs 3*(elnodes(2)-1)+1:3*elnodes(2)];
        E=elems(iel,3); I=elems(iel,4); A=elems(iel,5); Py=A*elems(iel,6); Mp=elems(iel,6)*elems(iel,7);
        Pi=FINT(1,iel);Mi=FINT(3,iel);Pj=FINT(4,iel);Mj=FINT(6,iel);
        [kt]=kelemp(nodexy,E,I,A,Py,Mp,Pi,Mi,Pj,Mj,tol);    %kt=globalized tangent stiffness matrix (Ziemian eq 8.5)
        Kt(eldofs,eldofs)=Kt(eldofs,eldofs) + kt;           %Assembly of Globalized Structure Stiffness Matrix
    end
    Keff=Kt(doffree,doffree);
    
    %-------------------Checking for System Limit State--------------------
    [V,D]=eig(Keff);                   %checking positive definiteness of stiffness matrix. Alternatively, one could check sign of v'*K11*v for some nonzero v.
    if find(diag(D<=0))                %Main idea is that at limit point, K11 ceases to be positive definite. Atleast one Eigenvalue <=0
        disp('Limit Point reached. Stiffness matrix is no longer positive definite. Program will terminate')
        break;
    elseif rcond(decomposition(Keff))<eps %checking ill-conditioning of stiffness matrix.
        disp('Stiffness Matrix has evolved into an Ill-conditioned Matrix. Collapse Mechanism Reached.Program will terminate')
        break;
    end
    
    du(doffree)=Keff\(dlamda*Pref(doffree)); %solving KD=F equation
    
    
    %--------------LOCAL STATE DETERMINATION---------------------------
    dfint=zeros(6,Nel);                %Initialization of localized internal force increment vector
    tau_array=[];                      %Initializing Tau array to store dlamda scale factor.
    dup=zeros(6,Nel);
    du_total=zeros(6,Nel);
    for iel=1:Nel
        elnodes=elems(iel,1:2);
        nodexy=nodes(elnodes,:);
        eldofs=3*(elnodes(1)-1)+1:3*elnodes(1);
        eldofs=[eldofs 3*(elnodes(2)-1)+1:3*elnodes(2)];
        E=elems(iel,3); I=elems(iel,4); A=elems(iel,5); Py=A*elems(iel,6); Mp=elems(iel,6)*elems(iel,7);
        Pi=FINT(1,iel);Mi=FINT(3,iel);Pj=FINT(4,iel);Mj=FINT(6,iel);
        [kt,T,G,ke,km]=kelemp(nodexy,E,I,A,Py,Mp,Pi,Mi,Pj,Mj,tol);
        ReqDef=T*du(eldofs);                    %global to local transformation
        
        pi=Pi/Py; mi=Mi/Mp; pj=Pj/Py; mj=Mj/Mp;
        phi_prev_i=pi^2 + mi^2 + 3.5*(pi*mi)^2; %yield condition state prior to load step:Near End
        phi_prev_j=pj^2 + mj^2 + 3.5*(pj*mj)^2; %yield condition state prior to load step:Far End
        
        %Checking for Elastic Unloading from Yield Surface/Storing
        %Plastic Deformations
        if abs(phi_prev_i-1)<=tol || abs(phi_prev_j-1)<=tol
            plasticmagnitude=(((G'*ke*G)^-1)*G'*ke)*ReqDef;
            if find(plasticmagnitude<0)
                disp('Attempted unload/reload. Program will terminate')
                unloading_ends=find(plasticmagnitude<0) %1==near end. 2==far end
                unloadingelem_no=iel
                return; %terminate program
            else
                dup(:,iel)=G*plasticmagnitude;   %storing element local axis plastic nodal displacements/rotations
            end
        end
        
        du_total(:,iel)=ReqDef;
        klocal=ke+km;
        dfint(:,iel)=klocal*ReqDef;              %localforce internal delta-force'increment' vector
        fint(:,iel)=FINT(:,iel) +dfint(:,iel);   %updating trial internal force vector
        
        %--UPDATING PHI ARRAY. SEE FORMAT OF PHI ARRAY IN SYSTEM INITIALIZATION----
        Pi_trial=fint(1,iel);Mi_trial=fint(3,iel);Pj_trial=fint(4,iel);Mj_trial=fint(6,iel);
        pi_trial=Pi_trial/Py; mi_trial=Mi_trial/Mp; pj_trial=Pj_trial/Py; mj_trial=Mj_trial/Mp;
        
        pm_trial_table=[pi_trial mi_trial pj_trial mj_trial];
        
        phi_trial_i=pi_trial^2 + mi_trial^2 + 3.5*(pi_trial*mi_trial)^2;                       %trial yield condition value for original surface:near end
        phi_trial_j=pj_trial^2 + mj_trial^2 + 3.5*(pj_trial*mj_trial)^2;                       %trial yield condition value for original surface:far end
        phi_tol_i=(1/1.01^2)*(pi_trial^2 + mi_trial^2 + (3.5/1.01^2)*(pi_trial*mi_trial)^2);   %trial yield cond.val.for expanded surface:near end
        phi_tol_j=(1/1.01^2)*(pj_trial^2 + mj_trial^2 + (3.5/1.01^2)*(pj_trial*mj_trial)^2);   %trial yield cond.val.for expanded surface:far end
        
        phi(:,iel)=[phi_trial_i;phi_prev_i;phi_tol_i;phi_trial_j;phi_prev_j;phi_tol_j];        %PHI ARRAY POPULATION
        
        phi_alarm_i=(1/0.99^2)*(pi_trial^2 + mi_trial^2 + (3.5/0.99^2)*(pi_trial*mi_trial)^2); %yield condition near end w.r.t. degenerated shrunk surface
        phi_alarm_j=(1/0.99^2)*(pj_trial^2 + mj_trial^2 + (3.5/0.99^2)*(pj_trial*mj_trial)^2); %yield condition far end w.r.t. degenerated shrunk surface
        
        if (phi_alarm_i <=(1+tol) && abs(phi_prev_i-1)<=tol) || ((phi_alarm_j <=(1+tol)) && abs(phi_prev_j-1)<=tol)
            %this means that trial pm point has gone under severe
            %degeneracy i.e. tangent vector to the original surface is
            %degenerate! analysis must terminate
            disp('Plastic flow vectors for some member ends might not be positioned along the surface properly. Internal force recovery might not be accurate.');
        end
        
        
        %----------LOAD-STEP SCALE FACTOR DETERMINATION------------------
        %dlamda scale factor determination for controlling surface drift/overshoot
        %populating tau array initialized earlier at the beginning of
        %state determination
        for j=1:2
            if phi((3*j)-2,iel)>(1+tol) && phi((3*j)-1,iel)<(1-tol)
                PM_prev=FINT([3*j-2 3*j],iel);
                PM_trial=fint([3*j-2 3*j],iel);
                A=elems(iel,5);Py=A*elems(iel,6);
                Mp=elems(iel,6)*elems(iel,7);
                [tau]=Elastic_Return(PM_prev,PM_trial,Py,Mp,tol); %Falsi Routine to scale dlamda for returning overshooting elastic points
                tau_array=[tau_array tau];
            elseif phi(3*j,iel)>(1+tol) && abs(phi((3*j)-1,iel)-1)<=tol
                PM_prev=FINT([3*j-2 3*j],iel);
                PM_trial=fint([3*j-2 3*j],iel);
                A=elems(iel,5);Py=A*elems(iel,6);
                Mp=elems(iel,6)*elems(iel,7);
                [tau]=Drift_Control(PM_prev,PM_trial,Py,Mp,tol); %Falsi Routine to scale dlamda for controlling surface drift to a defined tolerance
                tau_array=[tau_array tau];
            end
        end
    end
    
    %determining minimum scale factor (tau_govern) to scale incremental
    %force-displacement results
    if find(tau_array)
        tau_govern=min(tau_array(tau_array>0));
        dfint=tau_govern*dfint; %scaling current-step incremental forces
        du=tau_govern*du;
        dup=tau_govern*dup;
        du_total=tau_govern*du_total;
        dlamda=tau_govern*dlamda;
    end
    %Updating scaled total internal forces and return to surface
    for iel=1:Nel
        fint(:,iel)=FINT(:,iel) +dfint(:,iel); %here, fint is reevaluated after dfint has been scaled in close proximity to yield surface.
        Pi=FINT(1,iel);Mi=FINT(3,iel);Pj=FINT(4,iel);Mj=FINT(6,iel);
        A=elems(iel,5);Py=A*elems(iel,6); Mp=elems(iel,6)*elems(iel,7);
        pi=Pi/Py; mi=Mi/Mp; pj=Pj/Py; mj=Mj/Mp;
        phi_prev_i=pi^2 + mi^2 + 3.5*(pi*mi)^2;
        phi_prev_j=pj^2 + mj^2 + 3.5*(pj*mj)^2;
        Pi_trial=fint(1,iel);Mi_trial=fint(3,iel);Pj_trial=fint(4,iel);Mj_trial=fint(6,iel);
        pi_trial=Pi_trial/Py; mi_trial=Mi_trial/Mp; pj_trial=Pj_trial/Py; mj_trial=Mj_trial/Mp;
        
        pm_trial_table=[pi_trial mi_trial pj_trial mj_trial];
        
        phi_trial_i=pi_trial^2 + mi_trial^2 + 3.5*(pi_trial*mi_trial)^2;
        phi_trial_j=pj_trial^2 + mj_trial^2 + 3.5*(pj_trial*mj_trial)^2;
        phi_tol_i=(1/1.01^2)*(pi_trial^2 + mi_trial^2 + (3.5/1.01^2)*(pi_trial*mi_trial)^2);
        phi_tol_j=(1/1.01^2)*(pj_trial^2 + mj_trial^2 + (3.5/1.01^2)*(pj_trial*mj_trial)^2);
        
        phi(:,iel)=[phi_trial_i;phi_prev_i;phi_tol_i;phi_trial_j;phi_prev_j;phi_tol_j];  %Updating phi array
        
        %POINT PROJECTION OF SCALED LOCAL STATE ONTO YIELD SURFACE
        if point_projection_radial
            for x=1:2
                if abs(phi((3*x)-1,iel)-1)<=tol && phi((3*x),iel)<=(1+tol) && (phi((3*x)-2,iel)>(1+tol) || phi((3*x)-2,iel)<(1-tol))
                    if phi((3*x)-2,iel)>(1+tol)
                        %disp('drifted point above yield surface.')
                        flag=1;
                    elseif phi((3*x)-2,iel)<(1-tol)
                        %disp('drifted point below yield surface.')
                        flag=0;
                    end
                    pvalue=pm_trial_table(2*x -1);
                    mvalue=pm_trial_table(2*x);
                    [correc]=radial_projection(pvalue,mvalue,tol2,flag);
                    if flag==1
                        fint([3*x-2 3*x],iel)=correc*fint([3*x-2 3*x],iel);
                        pcorrec=correc*pvalue;
                        mcorrec=correc*mvalue;
                        phi([3*x-2 3*x],iel)= [pcorrec^2 + mcorrec^2 + 3.5*(pcorrec*mcorrec)^2;
                            (1/1.01^2)*(pcorrec^2 + mcorrec^2 + (3.5/1.01^2)*(pcorrec*mcorrec)^2)];
                        
                    else
                        fint([3*x-2 3*x],iel)=fint([3*x-2 3*x],iel)+[correc(1)*Py;correc(2)*Mp];
                        pcorrec=pvalue+correc(1);
                        mcorrec=mvalue+correc(2);
                        phi([3*x-2 3*x],iel)=[pcorrec^2 + mcorrec^2 + 3.5*(pcorrec*mcorrec)^2;
                            (1/1.01^2)*(pcorrec^2 + mcorrec^2 + (3.5/1.01^2)*(pcorrec*mcorrec)^2)];
                        
                    end
                end
            end
        elseif point_projection_normal
            for x=1:2
                if abs(phi((3*x)-1,iel)-1)<=tol && phi((3*x),iel)<=(1+tol) && (phi((3*x)-2,iel)>(1+tol) || phi((3*x)-2,iel)<(1-tol))
                    p1=pm_trial_table(2*x -1);
                    m1=pm_trial_table(2*x);
                    [msol,psol,alphasol]=normal_projection(m1,p1,Py,Mp);
                    fint([3*x-2 3*x],iel)=[psol*Py;msol*Mp];
                    phi([3*x-2 3*x],iel)=[msol^2 + psol^2 + 3.5*(msol*psol)^2;
                        (1/1.01^2)*(msol^2 + psol^2 + (3.5/1.01^2)*(msol*psol)^2)];
                end
            end
        end
        
        
        %identifying/recording first instance of a plastifying member end, node
        %number and lamda value.
        
        for j=1:2
            if abs(phi((3*j)-2,iel)-1)<=tol && phi((3*j)-1,iel)<(1-tol)
                node=elems(iel,j);
                hingetable=[hingetable;
                    iel node (lamda+dlamda)];
            end
        end
    end
    
    %-------UPDATING SYSTEM VARIABLES AFTER LOAD STEP------------------
    %updating system variables after satisfying yield surface
    %drift/overshoot
    FINT=fint;    %updating converged FINT array with current fint values
    FINT_history=[FINT_history;FINT];
    
    
    u=u+du;
    u_history=[u_history u];
    
    up=up+dup;
    up_history=[up_history;up];
    
    u_total=u_total +du_total;
    u_total_history=[u_total_history;u_total];
    
    lamda=lamda+dlamda;
    zz=zz+1;
    lamda_history(zz)=lamda;
    
    phi_history=[phi_history phi];
end



%--------------------------POST PROCESSING---------------------------------

%plotting node displacement/rotation (DSM- coordinates) vs Applied Load
%Ratio
freedof=4;   %DOF no to be specified by the Analyst for plotting purposes
xval=u_history(freedof,:);
yval=lamda_history;
figure
plot(xval,yval)
%xlim([0 10])
xlabel('in/rad')
ylabel('Load Ratio')
title(['DOF# ',num2str(freedof),' load vs displacement/rotation plot'])


%plotting Plastic as well as elastic member end displacement/rotation vs
%Applied Load Ratio
memberno=2; %To be specified.
endloc=2;   %To be specified. Near end=1, Far-end=2
dof=3;      %To be specified. For Axial disp, dof==1, For rotation, dof=3
xval1=u_total_history((3*endloc)-(3-dof):6:size(u_total_history(:,memberno)),memberno);
yval=lamda_history;
xval2=up_history((3*endloc)-(3-dof):6:size(up_history(:,memberno)),memberno);
figure
plot(abs(xval2),yval,'r',abs(xval1-xval2),yval,'*')
title('elastic-plastic member end disp/rotations vs load ratio Plot')
legend('u-plastic','u-elastic')
xlabel('in/rot')
ylabel('Applied Load Ratio')

%Plotting P-M point traces with P-M curve(Ziemian eq 10.18)
memberno=2;   %to be specified by the Analyst for plotting purposes
endloc=1; %To be specified. Near End==1, Far End==2.
if endloc==1
    location='near end';
elseif endloc==2
    location='far end';
end
Py=elems(memberno,5)*elems(memberno,6);
Mp=elems(memberno,6)*elems(memberno,7);
Ptrace=FINT_history((3*endloc)-2:6:size(FINT_history(:,memberno)),memberno);
Mtrace=FINT_history(3*endloc:6:size(FINT_history(:,memberno)),memberno);
%Ptrace_far=FINT_history(4:6:size(FINT_history(:,memberno)),memberno);
%Mtrace_far=FINT_history(6:6:size(FINT_history(:,memberno)),memberno);
M_curve=linspace(0,Mp,1000);
P_curve=Py*(sqrt((Mp^2-M_curve.^2)./(Mp^2+3.5.*M_curve.^2)));
%P_curve=Py*(sqrt((Mp^2-Mtrace_near.^2)./(Mp^2+3.5.*Mtrace_near.^2)));
figure
plot(abs(Mtrace),abs(Ptrace),'*',M_curve,P_curve)
xlabel('Moment/Kin')
ylabel('Axial Force/K')
title({'P-M trace history';['Member# ',num2str(memberno),'. location: ',location]})

%the script below filters out the hingetable to include first hinge
%formation events only.
%{
hingebook(1,:)=hingetable(1,:);
for i=2:size(hingetable(:,1))
      %if ~(find(hingebook(:,[1 2])==hingetable(i,1)) && find(hingebook(:,2)==hingetable(i,2)))
      if find(hingebook(:,1)==hingetable(i,1))
         if find(hingebook(
          
          hingebook(i,:)=hingetable(i,:);
      end
end
%}











