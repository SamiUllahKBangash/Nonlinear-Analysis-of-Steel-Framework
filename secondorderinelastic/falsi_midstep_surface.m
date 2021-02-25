function[tau_r] = falsi_midstep_surface(pval,mval,tol)
%Welcome to Falsi_midstep_surface. This function call takes trial secant
%step p and m coordinates in p-m space and computes isotropic expansion
%factor tau_r such that trial p,m values converge onto a yield surface
%corresponding to tau_r expansion factor. This fictitious yield surface
%will then enable us to compute secant G array for computing secant
%stiffness matrix later on. 

converge=false;
p=pval;
m=mval;

phi_trial_mid=p^2 + m^2 + 3.5*(p*m)^2;


if phi_trial_mid>(1+tol)
tau_a=0;
tau_b=0.01;
elseif phi_trial_mid<(1-tol)
tau_a=0;
tau_b=-0.01;
else
tau_r=0;
converge=true;

end

while ~converge
    
    phi_a=(1/(1 + tau_a)^2)*(p^2 + m^2 + (3.5/(1 + tau_a)^2)*(p^2*m^2));
    phi_b=(1/(1 + tau_b)^2)*(p^2 + m^2 + (3.5/(1 + tau_b)^2)*(p^2*m^2));
    
    tau_r= tau_b - ((phi_b -1)*((tau_a - tau_b)/(phi_a - phi_b)));
    
    phi_r=(1/(1 + tau_r)^2)*(p^2 + m^2 + (3.5/(1 + tau_r)^2)*(p^2*m^2));
    
    if abs(phi_r -1) <= tol
        converge=true;
        break;
    end
    
    if (phi_r -1)*(phi_b -1)>0
        tau_b=tau_r;
    elseif (phi_r-1)*(phi_a -1)>0
        tau_a=tau_r;
    end
end

end



    
    
    
    






