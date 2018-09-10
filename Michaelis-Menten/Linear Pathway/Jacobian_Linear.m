function m=Jacobian_Linear(x,Sin,Sout,kcat,Ks,Kp,Keq)
%Calculate Jacobian(d(dxdt)/dx) for a linear pathway with Michaelis-Menten
%kinetics
n=length(x);
d=zeros(n+1,2);
x_ext=[Sin;x(:);Sout];
for i=1:n+1
    d(i,:)=dMM(x_ext(i),x_ext(i+1),kcat(i),Ks(i),Kp(i),Keq(i))';
end
dvdx_s=sparse(2:n+1,1:n,d(2:n+1,1),n+1,n);
dvdx_p=sparse(1:n,1:n,d(1:n,2),n+1,n);
dvdx=dvdx_s+dvdx_p;
m=dvdx(1:n,:)-dvdx(2:n+1,:);
end