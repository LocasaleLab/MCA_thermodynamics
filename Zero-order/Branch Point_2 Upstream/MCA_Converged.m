function [fcc,J,dg]=MCA_Converged(k,K,Sin1,Sin2,Sout)
%Calculate flux control coefficients for a branched pathway
S_BP=((k(1)+k(2)-k(3))*K(1)*K(2)*K(3)*Sin1*Sin2...
    +(K(1)*K(2)*K(3)*Sin1*Sin2*(K(1)*K(2)*((k(1)+k(2)-k(3))^2)...
    *K(3)*Sin1*Sin2+4*k(3)*(K(1)*k(2)*Sin1+k(1)*K(2)*Sin2)*Sout))^0.5)/...
    (2*K(3)*(K(1)*k(2)*Sin1+k(1)*K(2)*Sin2));
if S_BP>0
    J1=k(1)-k(1)/K(1)*S_BP/Sin1;
    J2=k(2)-k(2)/K(2)*S_BP/Sin2;
    Mat=[1/J1 0 S_BP/(S_BP-K(1)*Sin1);...
        0 1/J2 S_BP/(S_BP-K(2)*Sin2);...
        1/(J1+J2) 1/(J1+J2) Sout/(K(3)*S_BP-Sout)];
    B=[1/J1 0 0;...
        0 1/J2 0;...
        1/(J1+J2) 1/(J1+J2) 0];
    J=[J1 J2 J1+J2];
    Q=[S_BP/Sin1 S_BP/Sin2 Sout/S_BP];
    dg=log(Q(:)./K(:))';
    fcc=(Mat'\B')';
else
    J=[-1 -1 -1];
    dg=[1 1 1];
    fcc=zeros(3,3);
end
end