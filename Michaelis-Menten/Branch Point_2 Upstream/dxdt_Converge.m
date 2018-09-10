function [dxdt,J]=dxdt_Converge(x,Sin1,Sin2,Sout,kcat,Ks,Kp,Keq)
%Calculate dxdt for converging pathway (2 upstream) with Michaelis-Menten kinetics
%x is concentration of the branch point metabolite
v=zeros(3,1); %rates of each reaction
v(1)=MM(Sin1,x,kcat(1),Ks(1),Kp(1),Keq(1));
v(2)=MM(Sin2,x,kcat(2),Ks(2),Kp(2),Keq(2));
v(3)=MM(x,Sout,kcat(3),Ks(3),Kp(3),Keq(3));
dxdt=v(1)+v(2)-v(3); %net changes of metabolites
J=Jacobian_Branch(x,Sin1,Sin2,Sout,kcat,Ks,Kp,Keq);
end