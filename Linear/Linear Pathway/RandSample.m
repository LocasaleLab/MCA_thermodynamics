nrxn=10;
nmet=nrxn-1;
nsamples=20000;
KSample=10.^normrnd(1,1,nrxn,nsamples);
kSample=10.^normrnd(1,1,nrxn,nsamples);
Sin=10.^normrnd(1,1,1,nsamples);
Sout=10.^normrnd(1,1,1,nsamples);

count=0;
for i=1:nsamples
    [S,f,fcc,dg]=MCA_Linear(kSample(:,i),KSample(:,i),Sin(i),Sout(i));
    if f>0
        count=count+1;
        SMat(count,:)=S(:)';
        fVec(count)=f;
        fccMat(count,:)=fcc(:)';
        dgMat(count,:)=dg(:)';
        kGood(count,:)=kSample(:,i)';
        KGood(count,:)=KSample(:,i)';
        giniVec(count)=GiniIndex(fcc);
    end
end
dgMax=max(dgMat');

figure;
for i=1:10
    subplot(2,5,i);
    dscatter(dgMax(:),fccMat(:,i));
    fccname=strcat('C^J_{v_{',num2str(i),'}}');
    title(fccname);
    xlabel('max\{g_i\}');
    ylabel(fccname);
    box on;
end
colormap('redbluecmap');
figure;
for i=1:10
    subplot(2,5,i);
    dscatter(dgMat(:,i),fccMat(:,i));
    fccname=strcat('C^J_{v_{',num2str(i),'}}');
    title(fccname);
    xlabel(strcat('g_{',num2str(i),'}'));
    ylabel(fccname);
    box on;
end
colormap('redbluecmap');