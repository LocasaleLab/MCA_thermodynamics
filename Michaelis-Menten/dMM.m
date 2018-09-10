function d=dMM(S,P,kcat,Ks,Kp,Keq)
%Calculate partial derivatives of reaction rate to substrate and product
%concentrations
%d(1): derivative against substrate
%d(2): derivative against product
d(1)=kcat*Kp*(Keq*Kp*Ks+Kp*P+Keq*Ks*P)/(Keq*(Kp*Ks+Ks*P+Kp*S)^2);
d(2)=-kcat*Kp*(Kp*Ks+Kp*S+Keq*Ks*S)/(Keq*(Kp*Ks+Ks*P+Kp*S)^2);
end