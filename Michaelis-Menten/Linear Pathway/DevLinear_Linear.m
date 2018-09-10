function [saturation,isinbound]=DevLinear_Linear(fcc,TDF,Sin,S,Sout,Ks,Kp)
%Quantify the deviation from linear kinetics
%saturation: maximal saturation term in all reactions
%isinbound: bool variable indicating if the FCCs satisfy constraints placed
%by the linear kinetics
isinbound=IsInBound_Linear(fcc,TDF);
n=length(fcc);
S_ext=[Sin;S(:);Sout];
sat_v=zeros(n,1);
for i=1:n
    sat_v(i)=Saturation(S_ext(i),S_ext(i+1),Ks(i),Kp(i));
end
saturation=max(sat_v);
end