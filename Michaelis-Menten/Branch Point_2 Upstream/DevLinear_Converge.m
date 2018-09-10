function [saturation,isinbound]=DevLinear_Converge(fcc,TDF,Sin1,Sin2,S,Sout,Ks,Kp)
%Quantify the deviation from linear kinetics
%saturation: maximal saturation term in all reactions
%isinbound: bool variable indicating if the FCCs satisfy constraints placed
%by the linear kinetics
isinbound=IsInBound_Converge(fcc,TDF);
sat_v=[Saturation(Sin1,S,Ks(1),Kp(1));Saturation(Sin2,S,Ks(2),Kp(2));...
    Saturation(S,Sout,Ks(3),Kp(3))];
saturation=max(sat_v);
end