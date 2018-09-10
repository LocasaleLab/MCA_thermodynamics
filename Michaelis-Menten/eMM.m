function e=eMM(S,P,kcat,Ks,Kp,Keq)
%Calculate normalized elasiticity coefficients
%e(1): elasticity coefficient of substrate
%e(2): elasticity coefficient of product
e(1)=-S/(P-Keq*S)*(Keq*Kp*Ks+Kp*P+Keq*Ks*P)/(Kp*Ks+Ks*P+Kp*S);
e(2)=P/(P-Keq*S)*(Kp*Ks+Kp*S+Keq*Ks*S)/(Kp*Ks+Ks*P+Kp*S);
end