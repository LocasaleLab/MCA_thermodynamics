function [saturation,isinbound]=DevLinear_Branch(fcc,TDF,Sin,S,Sout1,Sout2,Ks,Kp)
%Quantify the deviation from linear kinetics
%saturation: maximal saturation term in all reactions
%isinbound: bool variable indicating if the FCCs satisfy constraints placed
%by the linear kinetics
isinbound=IsInBound_Branch(fcc,TDF);
sat_v=[Saturation(Sin,S,Ks(1),Kp(1));Saturation(S,Sout1,Ks(2),Kp(2));...
    Saturation(S,Sout2,Ks(3),Kp(3))];
saturation=max(sat_v);
end