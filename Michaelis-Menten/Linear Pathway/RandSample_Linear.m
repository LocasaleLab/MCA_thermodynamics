nrxn=10;
nmet=nrxn-1;
nsamples=50000;
Keq=10.^normrnd(1,1,nrxn,nsamples);
kcat=10.^normrnd(1,1,nrxn,nsamples);
Ks=10.^normrnd(1,1,nrxn,nsamples);
Kp=10.^normrnd(1,1,nrxn,nsamples);
Sin=10.^normrnd(0,1,1,nsamples);
Sout=10.^normrnd(0,1,1,nsamples);

SMat=zeros(nsamples,nrxn-1);fVec=-ones(nsamples,1);fccMat=zeros(nsamples,nrxn);
dgMat=zeros(nsamples,nrxn);satV=zeros(nsamples,1);dev2boundV=zeros(nsamples,1);
count=0;
for i=1:nsamples
    dG_overall=log(Sout(i)/Sin(i)/prod(Keq(:,i)));
    if dG_overall>0
        b=Sin(i);
        Sin(i)=Sout(i);
        Sout(i)=b;
        kcat(:,i)=kcat(nrxn:-1:1,i);
        Ks(:,i)=Ks(nrxn:-1:1,i);
        Kp(:,i)=Kp(nrxn:-1:1,i);
        Keq(:,i)=Keq(nrxn:-1:1,i);
    end
    [S,f,fcc,dg]=SS_Linear(Sin(i),Sout(i),kcat(:,i),Ks(:,i),Kp(:,i),Keq(:,i));
    if f>0
        count=count+1;
        SMat(count,:)=S(:)';
        fVec(count)=f;
        fccMat(count,:)=fcc(:)';
        dgMat(count,:)=dg(:)'; 
        [s,b]=DevLinear_Linear(fcc(:),max(dg(:)),...
            Sin(i),S,Sout(i),Ks(:,i),Kp(:,i));
        satV(count)=s;
        dev2boundV(count)=b;
    end
end
dgMax=max(dgMat');
GP=find(fVec>0);
SMat=SMat(GP,:);fVec=fVec(GP,:);fccMat=fccMat(GP,:);dgMat=dgMat(GP,:);
satV=satV(GP);dev2boundV=dev2boundV(GP);dgMax=dgMax(GP);

figure;
[a,idx]=sort(satV,'descend');
for i=1:10
    subplot(2,5,i);
    scatter(dgMax(idx)',fccMat(idx,i),[],satV(idx),'Marker','.');
    fccname=strcat('C^J_{v_{',num2str(i),'}}');
    title(fccname);
    xlabel('max\{g_i\}');
    ylabel(fccname);
    box on;
end


figure;
scatter(satV,dev2boundV,'Marker','.')
xlabel('Maximal saturation term');
ylabel('Deviation from the bounds');
    
%categorize the saturation term into 10 bins and calculate the probability
%of satisfying the constraints in each bin
p_inbound_bin=zeros(10,1);
for i=1:10
    d2b_bin=dev2boundV(satV<=0.1*i & satV>=0.1*(i-1));
    p_inbound_bin(i)=sum(d2b_bin==0)/length(d2b_bin);
end