function [y,J,FCC,dG]=SS_Branch(Sin,Sout1,Sout2,kcat,Ks,Kp,Keq)
%Calculate steady state metabolite concentrations based on input/output
%substrate concentrations and enzyme parameters
%Calculate flux control coefficients (FCCs), flux (J) and reaction free energies
%based on the steady state metabolite concentrations
%-----------------------------------------------------------------------
%Calculate steady state concentrations
y0=mean([Sin;Sout1;Sout2]);
options = optimoptions('lsqnonlin','Display','off');
%y=fsolve(@(x) dxdt_Linear(x,Sin,Sout,kcat,Ks,Kp,Keq),y0,options);
[y,resnorm] = lsqnonlin(@(x) dxdt_Branch(x,Sin,Sout1,Sout2,kcat,Ks,Kp,Keq),y0,0,[],options);
J=zeros(3,1);
dG=zeros(3,1);
if resnorm>1e-6
    FCC=zeros(3,3);
    dG=zeros(3,1);
    J=-1*ones(3,1);
else
    %-----------------------------------------------------------------------
    %Calculate flux
    J(1)=MM(Sin,y,kcat(1),Ks(1),Kp(1),Keq(1));
    J(2)=MM(y,Sout1,kcat(2),Ks(2),Kp(2),Keq(2));
    J(3)=MM(y,Sout2,kcat(3),Ks(3),Kp(3),Keq(3));
    %-----------------------------------------------------------------------
    %Calculate FCCs
    e1=eMM(Sin,y,kcat(1),Ks(1),Kp(1),Keq(1));
    e2=eMM(y,Sout1,kcat(2),Ks(2),Kp(2),Keq(2));
    e3=eMM(y,Sout2,kcat(3),Ks(3),Kp(3),Keq(3));
    e_mat=[e1(2);e2(1);e3(1)];
    K=[1/J(1) 1/J(1);1/J(2) 0;0 1/J(3)];
    A=[K e_mat];
    B=[K zeros(3,1)];
    FCC=(A'\B')';
    %-----------------------------------------------------------------------
    %Calculate free energies
    dG(1)=log(y/Sin/Keq(1));
    dG(2)=log(Sout1/y/Keq(2));
    dG(3)=log(Sout2/y/Keq(3));
end
