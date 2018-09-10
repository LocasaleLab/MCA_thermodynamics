function m=Jacobian_Branch(x,Sin,Sout1,Sout2,kcat,Ks,Kp,Keq)
%Calculate Jacobian(d(dxdt)/dx) for a branch pathway with Michaelis-Menten
%kinetics
d1=dMM(Sin,x,kcat(1),Ks(1),Kp(1),Keq(1))';
d2=dMM(x,Sout1,kcat(2),Ks(2),Kp(2),Keq(2))';
d3=dMM(x,Sout2,kcat(3),Ks(3),Kp(3),Keq(3))';
m=d1(2)-d2(1)-d3(1);
end