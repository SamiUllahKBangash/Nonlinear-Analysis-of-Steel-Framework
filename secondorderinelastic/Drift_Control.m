function [tau]=Drift_Control(PM_prev,PM_trial,Py,Mp,tol)
p=PM_prev(1)/Py;
m=PM_prev(2)/Mp;
dp=(PM_trial(1)/Py)-p;
dm=(PM_trial(2)/Mp)-m;
tau_a=0;
tau_b=1;

converge=false;
j=1;
while ~converge %|| j==50
    phi_a=(1/1.01^2)*((p+tau_a*dp)^2 + (m+tau_a*dm)^2 + 3.5/(1.01^2)*((p+tau_a*dp)^2*(m+tau_a*dm)^2));
    phi_b=(1/1.01^2)*((p+tau_b*dp)^2 + (m+tau_b*dm)^2 + 3.5/(1.01^2)*((p+tau_b*dp)^2*(m+tau_b*dm)^2));
    
    tau_r=tau_b - ((phi_b -1)*(tau_a-tau_b))/(phi_a - phi_b);
    
    phi_r=(1/1.01^2)*((p+tau_r*dp)^2 + (m+tau_r*dm)^2 + 3.5/(1.01^2)*((p+tau_r*dp)^2*(m+tau_r*dm)^2));
    
    if abs(phi_r-1)<=tol
        converge=true;
        break;
    end
    if (phi_r -1)*(phi_b-1)>0
        tau_b=tau_r;
    elseif (phi_r -1)*(phi_a -1)>0
        tau_a=tau_r;
    end
    j=j+1;
end
tau=tau_r;
end
