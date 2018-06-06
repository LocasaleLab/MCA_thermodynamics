function [fcc,J,dg]=MCA_Branch(k,K,Sin,Sout1,Sout2)
%Calculate flux control coefficients for a branched pathway
S_BP=(Sin*K(1)*K(2)*K(3)*(k(1)-k(2)-k(3))+...
    ((K(1)*K(2)*K(3)*Sin*(k(1)-k(2)-k(3)))^2+...
    4*k(1)*K(2)*K(3)*(K(1)*k(2)*K(3)*Sin*Sout1+...
    K(1)*K(2)*k(3)*Sin*Sout2))^0.5)/(2*k(1)*K(2)*K(3));
J1=k(2)-k(2)/K(2)*Sout1/S_BP;
J2=k(3)-k(3)/K(3)*Sout2/S_BP;
Mat=[1/(J1+J2) 1/(J1+J2) S_BP/(S_BP-K(1)*Sin);...
    1/J1 0 Sout1/(K(2)*S_BP-Sout1);...
    0 1/J2 Sout2/(K(3)*S_BP-Sout2)];
B=[1/(J1+J2) 1/(J1+J2) 0;1/J1 0 0;0 1/J2 0];
J=[J1+J2 J1 J2];
Q=[S_BP/Sin Sout1/S_BP Sout2/S_BP];
dg=log(Q(:)./K(:))';
fcc=(Mat'\B')';
end