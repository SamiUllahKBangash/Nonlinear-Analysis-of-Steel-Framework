%the routines developed below help illustrate the performance of Euler,
%RK-2, Newton Raphson Method (NRM) and Modified Newton Raphson Method
%(MNRM) in running a simplified second order elastic analysis of a
%reference 1 bay frame.

%the script requires external computation of condensed scalar elastic
%stiffness constant for the DOF under consideration. L depends on the frame
%geometry and c is assumed to be given
clc;
clear all;
ke=41.62; %The value of ke needs to be defined from prior first order elastic analysis of the given frame.
L=180;    %the height of the single storey frame
c=0.025;   %proportional loading constant H=c*P
lamdacr=c*ke*L;  %critical lateral load defined as point of zero secant stiffness
i=1;             
disp(i)=0;       %array to store top node displacement_euler
H(i)=0;          %array to store top node laterally applied load_euler
n=5;           %no of load steps defined by user
dlamda=lamdacr/n;   %load increments assumed as uniformly distributed

for i=2:n+1
    ddisp=(((ke*((c*L)/(c*L + disp(i-1)))^2))^-1)*dlamda;
    disp(i)=disp(i-1) + ddisp;
    H(i)=H(i-1) + dlamda;
end





%the code below computes runge kutta second order midpoint version of
%simplified second order analysis of the same frame for comparison with
%euler method above
ke=41.62; %The value of ke needs to be defined from prior first order elastic analysis of the given frame.
L=180;    %the height of the single storey frame
c=0.025;   %proportional loading constant H=c*P
lamdacr=c*ke*L;  %critical lateral load defined as the point of zero secant stiffness
i=1;             
disp2(i)=0;       %array to store top node displacement_RK-2
H2(i)=0;          %array to store top node laterally applied loa_RK-2
%n=100;           %RK-2 steps
dlamda=lamdacr/n;   %load increments assumed as uniformly distributed

for i=2:n+1
    ddisp2_predictor= (((ke*((c*L)/(c*L + disp2(i-1)))^2))^-1)*dlamda/2;
    disp2(i)=disp2(i-1) + ddisp2_predictor;
    K2bar=ke*((c*L)/(c*L + disp2(i)))^2;
    ddisp2_corrector=(K2bar^-1)*dlamda;
    disp2(i)=disp2(i-1) + ddisp2_corrector;
    H2(i)=H(i-1) + dlamda;
end






%the code below executes standard/modified NRM  for the above frame
%problem and compares results with euler and RK-2
tic
ke=41.62;                  %The value of ke needs to be defined from prior first order elastic analysis of the given frame.
L=180;                     %the height of the single storey frame
c=0.025;                   %proportional loading constant H=c*P
lamdacr=c*ke*L;            %critical lateral load defined as point of zero secant stiffness
i=1;             
disp3(i,i)=0;              %array to store top node displacement
D(i)=0;                    %converged displacement array
H3(i)=0;                   %array to store top node laterally applied load
Fint2(i,i)=0;              %array to store internal force for iterative history
%n=5;                      %(M)/NRM steps
dlamda_1=0.97*lamdacr/n;   %load increments assumed as uniformly distributed (reduced from full lamdacr for stable convergence)
tol=10^-6;                 %tolerance parameter. Residual convergence is affected by the magnitude
lamda=0;                   %loading variable
Href=1;                    %reference load
flagMNRM=false;            %flag to activate modified NRM method
Res(i,i)=0;                %variable to store force imbalance/residuals for each substep iteration

for i=2:n+1                %loop for loading increments in n steps
 j=1;                      %iteration counter
 converge=false;           %flag to check for convergence. set to false at start of each load step
 dtotal=0;                 %variable to store iterative displacements in each load step
 while ~converge
     if j==1
         dlamda=dlamda_1;  %loading increment at start of iteration for each load step 
         res_1=0;          %residual at start of iteration assumed zero for each load step
         d=D(i-1);         %converged displacement from previous load step
         kt_MNRM=(ke*((c*L)/(c*L + d))^2);  %modified NRM tangent formed at converged displacement state from previous load step
     else
         dlamda=0;         %loading increment for subsequent iterations within the load step, = 0 for (M)/NRM method in subsequent iterations
         d=disp3(j-1,i);   %displacement state from previous iteration with in load step
         
     end
          if flagMNRM      %check for activating MNRM tangent: subject to user discretion
              kt=kt_MNRM;
          else
              kt=(ke*((c*L)/(c*L + d))^2);  %standard NRM tangent
          end
     ddbar=(kt^-1)*Href;                  %See Ziemian eq 12.14a
     ddbarbar=(kt^-1)*res_1;              %See Ziemian eq 12.14b
     ddtotal=(dlamda*ddbar) + ddbarbar;   %See Ziemian eq 12.15
     dtotal=dtotal+ddtotal;               %See Ziemian eq 12.10
     disp3(j,i)=D(i-1)+dtotal;            %See Ziemian eq 12.10
     Fint2(j,i)=(ke*disp3(j,i))/(1+ (disp3(j,i)/(c*L)));  %array to store internal force for each iteration with in load step
     lamda=lamda+dlamda;                  %external force variable with loading increments at each iteration
     res_1=lamda - Fint2(j,i);            %Force imbalance/Residual
     Res(j,i)=res_1;                      %array to store Residuals at each load step iteration for investigation purposes  
         if norm(Res(j,i))<=tol || j==100  %check for convergence based on tolerance param. Or whether max allowable iterations exhausted or not
             D(i)=disp3(j,i);                       %variable to store converged displacement against each load step
             H3(i)=lamda;                  %variable to store converged load step.
             converge=true;                %boolean to exit iterations and begin new load step
         end
     j=j+1;                                %iteration counter with in load step
 end
end
toc

%define exact force displacement curve here
x=linspace(0,145.5,1000);
y=(ke*x)./(1+(x./(c*L)));



%POST PROCESSING/Visualization

p=1;
for z=2:n+1     %this script is written to trace the entire characteristic zig-zag iteration path for (M)/NRM method to help in visualization
    for j=1:find(disp3(:,z)~=0,1,'last')
         y2(p)=H3(z);
         x2(p)=disp3(j,z);
         y2(p+1)=Fint2(j,z);
         x2(p+1)=disp3(j,z);
         p=p+2;
    end
end
figure
plot(disp,H,'r',disp2,H2,'g',D,H3,'*',x,y,'b')
xlabel('lateral displacement/in')
ylabel('lateral load at top node/kips')
title('euler vs RK-2 vs NRM: simplified second order analysis')
legend('euler','RK-2','NRM','exact')

figure
plot(x,y,'b',D,H3,'*',[0 x2],[0 y2],'r')
xlabel('lateral displacement/in')
ylabel('lateral load at top node/Kips')
title('(M)/NRM iteration visualization')
legend('exact equilibrium curve','converged (M)/NRM solution','iteration history')


