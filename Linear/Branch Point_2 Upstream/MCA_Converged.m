function [fcc,J,dg]=MCA_Converged(k,K,Sin1,Sin2,Sout)
%Calculate flux control coefficients for a branched pathway
S_BP=(k(1)*Sin1+k(2)*Sin2+k(3)/K(3)*Sout)/(k(2)/K(2)+k(3)+k(1)/K(1));
J1=k(1)*Sin1-k(1)/K(1)*S_BP;
J2=k(2)*Sin2-k(2)/K(2)*S_BP;
Mat=[1/J1 0 -k(1)/K(1)*S_BP/J1;0 1/J2 -k(2)/K(2)*S_BP/J2;1/(J1+J2) 1/(J1+J2) k(3)*S_BP/(J1+J2)];
B=[1/J1 0 0;0 1/J2 0;1/(J1+J2) 1/(J1+J2) 0];
J=[J1+J2 J1 J2];
Q=[S_BP/Sin1 S_BP/Sin2 Sout/S_BP];
dg=log(Q(:)./K(:))';
fcc=(Mat'\B')';
end