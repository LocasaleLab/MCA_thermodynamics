function [fcc,J,dg]=MCA_Branch(k,K,Sin,Sout1,Sout2)
%Calculate flux control coefficients for a branched pathway
S_BP=(k(1)*Sin+k(2)/K(2)*Sout1+k(3)/K(3)*Sout2)/(k(2)+k(3)+k(1)/K(1));
J1=k(2)*S_BP-k(2)/K(2)*Sout1;
J2=k(3)*S_BP-k(3)/K(3)*Sout2;
Mat=[1/(J1+J2) 1/(J1+J2) -k(1)/K(1)/(J1+J2);1/J1 0 k(2)/J1;0 1/J2 k(3)/J2];
B=[1/(J1+J2) 1/(J1+J2) 0;1/J1 0 0;0 1/J2 0];
J=[J1+J2 J1 J2];
Q=[S_BP/Sin Sout1/S_BP Sout2/S_BP];
dg=log(Q(:)./K(:))';
fcc=(Mat'\B')';
end