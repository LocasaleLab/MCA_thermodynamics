function s=Saturation(S,P,Ks,Kp)
%Calculation the saturation term of a Michaelis-Menten reaction
s=1-1/(1+S/Ks+P/Kp);
%s=S/Ks/(1+S/Ks+P/Kp);
end