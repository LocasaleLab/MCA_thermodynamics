function v=MM(S,P,kcat,Ks,Kp,Keq)
%Kinetic law of Michaelis-Menten mechanism
%Parameters
%S: concentration of substrate
%P: concentration of product
%kcat,Ks,Kp,Keq: catalysis constant,Km for substrate and product,
%equilibrium constant
v=kcat*S/Ks/(1+S/Ks+P/Kp)*(1-P/S/Keq);