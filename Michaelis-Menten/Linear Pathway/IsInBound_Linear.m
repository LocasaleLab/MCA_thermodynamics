function a=IsInBound_Linear(fcc,TDF)
%Justify if FCCs satisfy the constraints in the same model with linear
%kinetics
n=length(fcc);
et=exp(TDF);
x=-fcc;
x(1)=fcc(1);
for i=1:n
    x_lb(i)=-(et^(i-1))/(1-et^n);
end
x_lb(1)=1-et;
delta_x=x(:)-x_lb(:);
delta_x(delta_x>0)=0;
%a=delta_x;
a=norm(delta_x);
end