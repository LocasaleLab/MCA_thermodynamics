function [dxdt,J]=dxdt_Linear(x,Sin,Sout,kcat,Ks,Kp,Keq)
%Calculate dxdt for linear pathway with Michaelis-Menten kinetics
n=length(x);
v=zeros(n+1,1); %rates of each reaction
v(1)=MM(Sin,x(1),kcat(1),Ks(1),Kp(1),Keq(1));
for i=2:n
    v(i)=MM(x(i-1),x(i),kcat(i),Ks(i),Kp(i),Keq(i));
end
v(n+1)=MM(x(n),Sout,kcat(n+1),Ks(n+1),Kp(n+1),Keq(n+1));
dxdt=v(1:n)-v(2:n+1); %net changes of metabolites
J=Jacobian_Linear(x,Sin,Sout,kcat,Ks,Kp,Keq);
end
