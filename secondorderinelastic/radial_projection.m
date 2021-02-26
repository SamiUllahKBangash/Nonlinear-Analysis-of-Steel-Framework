function [correc]=radial_projection(pvalue,mvalue,tol2,flag)
%------------------------option b------------------------------------------
if flag==1
    p=0;
    m=0;
    dp=pvalue;
    dm=mvalue;
    tau_a=0;
    tau_b=1;
else
    
    p=pvalue;
    m=mvalue;
    dp=sqrt(pvalue^2 + mvalue^2)*pvalue;
    dm=sqrt(pvalue^2 + mvalue^2)*mvalue;
    tau_a=0;
    tau_b=1;
end

%tol=0.01; %tolerance parameter for convergence to phi==1
converge=false;
j=1;
while ~converge %|| j==50
    phi_a=(p+tau_a*dp)^2 + (m+tau_a*dm)^2 + 3.5*(p+tau_a*dp)^2*(m+tau_a*dm)^2;
    phi_b=(p+tau_b*dp)^2 + (m+tau_b*dm)^2 + 3.5*(p+tau_b*dp)^2*(m+tau_b*dm)^2;
    
    tau_r=tau_b - ((phi_b -1)*(tau_a-tau_b))/(phi_a - phi_b);
    
    phi_r=(p+tau_r*dp)^2 + (m+tau_r*dm)^2 + 3.5*(p+tau_r*dp)^2*(m+tau_r*dm)^2;
    
    if (abs(phi_r-1))<=tol2
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

if flag==1
    correc=tau_r;
else
    correc=tau_r*[dp;dm];
end
end
