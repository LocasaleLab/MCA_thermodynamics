function [dxdt,J]=dxdt_Branch(x,Sin,Sout1,Sout2,kcat,Ks,Kp,Keq)
%Calculate dxdt for branch pathway (2 downstream) with Michaelis-Menten kinetics
%x is concentration of the branch point metabolite
v=zeros(3,1); %rates of each reaction
v(1)=MM(Sin,x,kcat(1),Ks(1),Kp(1),Keq(1));
v(2)=MM(x,Sout1,kcat(2),Ks(2),Kp(2),Keq(2));
v(3)=MM(x,Sout2,kcat(3),Ks(3),Kp(3),Keq(3));
dxdt=v(1)-v(2)-v(3); %net changes of metabolites
J=Jacobian_Branch(x,Sin,Sout1,Sout2,kcat,Ks,Kp,Keq);
end