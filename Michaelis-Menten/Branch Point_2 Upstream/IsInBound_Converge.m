function a=IsInBound_Converge(fcc,TDF)
%Justify if FCCs satisfy the constraints in the same model with linear
%kinetics
fcc_matrix=reshape(fcc,3,3);
et=exp(TDF);
x=[fcc_matrix(1,1);fcc_matrix(1,2);-fcc_matrix(1,3);fcc_matrix(2,1);...
    fcc_matrix(2,2);-fcc_matrix(2,3);-fcc_matrix(3,3)];
x_lb=[1-et;-et/(1-et);-et/(1-et);-et/(1-et);1-et;-et/(1-et);-et];
delta_x=x-x_lb;
delta_x(delta_x>0)=0;
a=norm(delta_x);
end