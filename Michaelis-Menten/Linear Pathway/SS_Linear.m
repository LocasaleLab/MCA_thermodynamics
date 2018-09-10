function [y,J,FCC,dG]=SS_Linear(Sin,Sout,kcat,Ks,Kp,Keq)
%Calculate steady state metabolite concentrations based on input/output
%substrate concentrations and enzyme parameters
%Calculate flux control coefficients (FCCs), flux (J) and reaction free energies
%based on the steady state metabolite concentrations
%-----------------------------------------------------------------------
%Calculate steady state concentrations
n=length(kcat)-1;
y0=((Sin*Sout)^0.5)*ones(n,1);
options = optimoptions('lsqnonlin','Display','off');
%y=fsolve(@(x) dxdt_Linear(x,Sin,Sout,kcat,Ks,Kp,Keq),y0,options);
[y,resnorm] = lsqnonlin(@(x) dxdt_Linear(x,Sin,Sout,kcat,Ks,Kp,Keq),y0,zeros(n,1),[],options);
if resnorm>1e-6
    FCC=zeros(n+1,1);
    dG=zeros(n+1,1);
    J=-1;
else
    %-----------------------------------------------------------------------
    %Calculate FCCs
    x_ext=[Sin;y(:);Sout]; %vector storing concentrations of all metabolites (Sin,Sout and intermediates)
    e=zeros(n+1,2);
    for i=1:n+1
        e(i,:)=eMM(x_ext(i),x_ext(i+1),kcat(i),Ks(i),Kp(i),Keq(i));
    end
    idx_i=[1:n+1;1:n+1];
    idx_i=idx_i(:);
    idx_j=[1:n;1:n];
    idx_j=[0;idx_j(:);0];
    e_vec=e';
    e_vec=e_vec(:);
    e_mat=sparse(idx_i(idx_j>0),idx_j(idx_j>0),e_vec(idx_j>0));
    A=[e_mat';ones(1,n+1)];
    b=[zeros(n,1);1];
    FCC=A\b;
    %-----------------------------------------------------------------------
    %Calculate free energies
    dG=log(x_ext(2:n+2)./x_ext(1:n+1)./Keq);
    %-----------------------------------------------------------------------
    %Calculate flux
    JVec=zeros(n+1,1);
    S_ext=[Sin;y(:);Sout];
    for i=1:n+1
        JVec(i)=MM(S_ext(i),S_ext(i+1),kcat(i),Ks(i),Kp(i),Keq(i));
    end
    J=mean(JVec);
    if (max(JVec)-min(JVec))/J>1e-4
        J=-1;
    end
end
