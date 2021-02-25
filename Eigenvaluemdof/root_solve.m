x=linspace(0,40,100);
y=x.*cos(x)- sin(x);

m=1;
for i=1:99
    
    if y(i+1)/y(i) <0
        sol(m)=x(i);
        m=m+1;
    end
end

sol
sol2=pi*(sol.^-1)
plot(x,y)
