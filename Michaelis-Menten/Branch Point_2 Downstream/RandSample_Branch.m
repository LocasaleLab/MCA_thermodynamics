nsamples=50000;
Keq=10.^normrnd(1,1,3,nsamples);
kcat=10.^normrnd(1,1,3,nsamples);
Ks=10.^normrnd(1,1,3,nsamples);
Kp=10.^normrnd(1,1,3,nsamples);
Sin=10.^normrnd(0,1,1,nsamples);
Sout1=10.^normrnd(0,1,1,nsamples);
Sout2=10.^normrnd(0,1,1,nsamples);

count=0;
SMat=zeros(nsamples,1);fMat=-ones(nsamples,3);fccMat=zeros(nsamples,9);
dgMat=zeros(nsamples,3);satV=zeros(nsamples,1);dev2boundV=zeros(nsamples,1);
for i=1:nsamples
    [S,f,fcc,dg]=SS_Branch(Sin(i),Sout1(i),Sout2(i),kcat(:,i),Ks(:,i),Kp(:,i),Keq(:,i));
    if min(f)>0
        count=count+1;
        SMat(count,:)=S(:)';
        fMat(count,:)=f(:)';
        fccMat(count,:)=fcc(:)';
        dgMat(count,:)=dg(:)';      
        [satV(count),dev2boundV(count)]=DevLinear_Branch(fcc(:),max(dg(:)),...
            Sin(i),S,Sout1(i),Sout2(i),Ks(:,i),Kp(:,i));
    end
end
dgMax=max(dgMat');
GP=find(abs((fMat*[1 -1 -1]')./min(fMat')')<1e-4 & min(fMat')'>0);
SMat=SMat(GP,:);fMat=fMat(GP,:);fccMat=fccMat(GP,:);dgMat=dgMat(GP,:);
satV=satV(GP);dev2boundV=dev2boundV(GP);dgMax=dgMax(GP);

for i=1:3
    for j=1:3
        subplot(3,3,(i-1)*3+j);
        %dscatter(dgMat(j,GP)',fccMat((j-1)*3+i,GP)'); %FCC(J_i,v_j)
        dscatter(dgMax(:),fccMat(:,(j-1)*3+i));
        ylim([-2 2]);
        title(strcat('C_{v_',num2str(j),'}^{J_',num2str(i),'}'));
        xlabel('max\{g_i\}');
        ylabel('FCC');
        box on;
    end
end
colormap('redbluecmap');

figure;
for i=1:3
    for j=1:3
        subplot(3,3,(i-1)*3+j);
        dscatter(dgMat(:,j),fccMat(:,(j-1)*3+i)); %FCC(J_i,v_j)
        ylim([-2 2]);
        title(strcat('C_{v_',num2str(j),'}^{J_',num2str(i),'}'));
        xlabel(strcat('g_',num2str(j)));
        ylabel('FCC');
        box on;
    end
end
colormap('redbluecmap');

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