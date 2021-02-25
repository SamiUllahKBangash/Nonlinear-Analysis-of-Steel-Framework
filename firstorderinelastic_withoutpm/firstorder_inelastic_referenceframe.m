


E=29e3;  %ksi
fy=36;   %ksi
rho=0.0001;   %hardening parameter
W=100;    %axial load kips
dh=1;       %step load
zz=1; %counter to update situation
Mint=zeros(2,4);
 disp_node2x=0;
 H=0;
 lateral_force(zz)=0;
 lateral_disp_node2x(zz)=0;

 for situation=1:1

%plastic flag variable. 

%the limitation of this code is that it assumes you already know the
%sequence of plastic hinge formation i.e. which hinges will plastify first.
%This is evident from the situations described below. Also, the code runs
%only for the given frame geometry. The code needs to be adapted to Prakash
%code for making it more general.
if situation==1
        elem1ith=0;
        elem1jth=0;
        elem2ith=0;
        elem2jth=0;
        elem3ith=0;
        elem4ith=0;
elseif situation==2
        elem1ith=0;
        elem1jth=0;
        elem2ith=0;
        elem2jth=0;
        elem3ith=0;
        elem4ith=1;
elseif situation==3
        elem1ith=0;
        elem1jth=1;
        elem2ith=0;
        elem2jth=0;
        elem3ith=0;
        elem4ith=1;
elseif situation==4
        elem1ith=0;
        elem1jth=1;
        elem2ith=0;
        elem2jth=0;
        elem3ith=1;
        elem4ith=1;
end

%Column details W14x132
lc=120;   %in
ac=38.8;  %in^2
ic=1530;  %in^4
zc=234;   %in^3
mpc=fy*zc;  %kin

%beam details W18x76
lb=120; %in
ib=1330; %in^4
ab=22.3;  %in^2
zb=163;  %in^3
mpb=fy*zb;   %kin


%element stiffness matrices


%beam stiffness (unplastified)
kb_1=rho*E*ib/lb*[ ab/ib 0 0;     %perfectly elastic portion
                    0    4 2;
                    0    2 4];
kb_2=(1-rho)*E*ib/lb*[ ab/ib 0 0;  %EPP portion (subject to change)
                        0    4 2;
                        0    2 4];
Kb1=kb_1+kb_2;                      %combined beam stiffness

%beam stiffness (plastified)
kb_1=rho*E*ib/lb*[ ab/ib 0 0;     %perfectly elastic portion
                    0    4 2;
                    0    2 4];
kb_2=(1-rho)*E*ib/lb*[ ab/ib 0 0;  %EPP portion (subject to change)
                        0    0 0;
                        0    0 3];
Kb2=kb_1+kb_2;                      %combined beam stiffness

if elem4ith==0
   Kb_elem_4=Kb1;
else
     Kb_elem_4=Kb2;
end
if elem3ith==0
    Kb_elem_3=Kb1;
else
    Kb_elem_3=Kb2;
end


%column stiffness (unplastified)
kc_1=rho*E*ic/lc*[ ac/ic  0 0;     %perfectly elastic portion
                     0    4 2;
                     0    2 4];
kc_2=(1-rho)*E*ic/lc*[ ac/ic  0 0;   %EPP portion (subject to change)
                       0      4 2;
                       0      2 4];
Kc1=kc_1+kc_2;                        %combined column stiffness

%column stiffness (plastified ith end)

kc_1=rho*E*ic/lc*[ ac/ic  0 0;     %perfectly elastic portion
                     0    4 2;
                     0    2 4];
kc_2=(1-rho)*E*ic/lc*[ ac/ic  0 0;   %EPP portion (subject to change)
                       0      0 0;
                       0      0 3];
Kc2_ith=kc_1+kc_2; 

%column stiffness (plastified jth end)

kc_1=rho*E*ic/lc*[ ac/ic  0 0;     %perfectly elastic portion
                     0    4 2;
                     0    2 4];
kc_2=(1-rho)*E*ic/lc*[ ac/ic  0 0;   %EPP portion (subject to change)
                       0      3 0;
                       0      0 0];
Kc2_jth=kc_1+kc_2; 


%column stiffness (plastified ith and jth end)
kc_1=rho*E*ic/lc*[ ac/ic  0 0;     %perfectly elastic portion
                     0    4 2;
                     0    2 4];
kc_2=(1-rho)*E*ic/lc*[ ac/ic  0 0;   %EPP portion (subject to change)
                       0      0 0;
                       0      0 0];
Kc2_ij=kc_1+kc_2; 


if elem1ith==0 && elem1jth==0
    Kc_elem_1=Kc1;
elseif elem1ith==1 && elem1jth==0
    Kc_elem_1=Kc2_ith;
elseif elem1ith==0 && elem1jth==1
    Kc_elem_1=Kc2_jth;
elseif elem1ith==1 && elem1jth==1
    Kc_elem_1=Kc2_ij;
end
    
 if elem2ith==0 && elem2jth==0
    Kc_elem_2=Kc1;
elseif elem2ith==1 && elem2jth==0
    Kc_elem_2=Kc2_ith;
elseif elem2ith==0 && elem2jth==1
    Kc_elem_2=Kc2_jth;
elseif elem2ith==1 && elem1jth==1
    Kc_elem_2=Kc2_ij;   
 end

%transformation matrix

at=[1  0   0 -1  0    0;             %element transformation matrix
    0 1/lc 1 0 -1/lc 0;
    0 1/lc 0 0 -1/lc 1];
tr2=[0 1 0;                            
   -1 0 0;
    0 0 1];
Tr2=[tr2 zeros(3);                   %coordinate transformation matrix for 2nd story column
    zeros(3) tr2];

tr1=[0 -1 0;
   1 0 0;
    0 0 1];
Tr1=[tr1 zeros(3);                    %coordinate transformation matrix for 1st story column
    zeros(3) tr1];


%transformed K


Kglobe_col1=Tr1'*(at'*Kc_elem_1*at)*Tr1;   %globalised 1st story column stiffness
Kglobe_col2=Tr2'*(at'*Kc_elem_2*at)*Tr2;   %globalised 2nd story column stiffness
Kglobe_beam3=at'* Kb_elem_3*at;              %globalised beam stiffness
Kglobe_beam4=at'* Kb_elem_4*at;


Ktotal=zeros(15);                   %global stiffness    
eft=[1 2 3 13 14 15;      %elem1          %element dof table  row wise
       1   2  3 4 5 6;     %elem2          
       4   5  6 7 8 9;      %elem3 
       1   2  3 10 11 12];   %elem4 
   
%assembly algorithm 
for s=1:4
    EFT=eft(s,:);
    if s==1 
        klocal=Kglobe_col1;
    elseif s==2
        klocal=Kglobe_col2;
    elseif s==3
        klocal=Kglobe_beam3;
    elseif s==4
        klocal=Kglobe_beam4;
    end
   for i=1:6
    for j=1:6
       Ktotal(EFT(i),EFT(j))=Ktotal(EFT(i),EFT(j))+klocal(i,j);
    end
   end
end

 K_original=Ktotal;    %original stiffness _unmodified
 Ktotal_mod=Ktotal;    %modified constrained stiffness (non singular)
 cst=[8 11 13 14 15];
 
 %stiffness constraining algorithm
 for c=1:5
 Ktotal_mod(cst(c),:)=0*Ktotal(cst(c),:);
  Ktotal_mod(:,cst(c))=0*Ktotal(:,cst(c));
  Ktotal_mod(cst(c),cst(c))=1;
 end
 
 
 %DH=[0.3*dh;-dh;0;0.3*dh;-dh;0;0;0;0;0;0;0;0;0;0];   %global force vector RHS matrix (modified)
 DH=[dh;0;0;dh;0;0;0;0;0;0;0;0;0;0;0];
 dd=Ktotal_mod\DH;                       %displacement vector global   
 
 %state space determination
 elem_force=zeros(6,4);
 for s=1:4
     EFT=eft(s,:);
     if s==1 
      klocal=at'*Kc_elem_1*at;
      TR=Tr1;
     elseif s==2
         klocal=at'*Kc_elem_2*at;
         TR=Tr2;
     elseif s==3
         klocal=at'*Kb_elem_3*at;
         TR=eye(6);
     elseif s==4
         klocal=at'*Kb_elem_4*at;
         TR=eye(6);
     end
      for i=1:6
          delem(i)=dd(EFT(i));
      end
      elem_force(:,s)=klocal*TR*delem';
 end
 %element moments only
 elem_moment=zeros(2,4);
 for z=1:4
     elem_moment(1,z)=elem_force(3,z);
     elem_moment(2,z)=elem_force(6,z);
 end
 
 
 %compute scale factors for each member end

 

  sf=zeros(2,4);
 for i=1:4
     if i==1
         Mp=mpc;
     elseif i==2
         Mp=mpc;
     elseif i==3
         Mp=mpb;
     else
         Mp=mpb;
     end
  for j=1:2 
     
      %if Mint(j,i)/elem_moment(j,i) < 0
         %print('elastic unloading occuring!')
       %  return
     if abs(elem_moment(j,i))>0.001
              sf(j,i)=(Mp-abs(Mint(j,i)))/abs(elem_moment(j,i));
     end
  end
 end
 zz=zz+1;
 %scale factor that will be used to scale up the current reference load, nodal displacements,elem_moment
 SF_govern=min(sf(sf > 0));
 Mint=Mint+SF_govern*(elem_moment);
 disp_node2x=disp_node2x+SF_govern*dd(4);
 H=H+SF_govern*dh;
 lateral_force(zz)=H;
 lateral_disp_node2x(zz)=disp_node2x;
 end
 
plot(lateral_disp_node2x,0.3*lateral_force)
xlabel('node 2 horizontal displacement/in')
ylabel('horizontal lateral load/kips')
xlim([0 24])

