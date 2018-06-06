NSample=20000;
k=10.^normrnd(1,1,3,NSample);
K=10.^normrnd(1,1,3,NSample);
Sin=10.^normrnd(1,1,1,NSample);
Sout1=10.^normrnd(1,1,1,NSample);
Sout2=10.^normrnd(1,1,1,NSample);

for i=1:NSample
    [fcc,J,dg]=MCA_Branch(k(:,i),K(:,i),Sin(i),Sout1(i),Sout2(i));
    fccMat(:,i)=fcc(:);
    JMat(:,i)=J(:);
    dgMat(:,i)=dg(:);
end

dgMax=max(dgMat);
dgMin=min(dgMat);
GP=find(min(JMat)>0); %GoodPos, which means fluxes through each reaction is positive
for i=1:3
    for j=1:3
        subplot(3,3,(i-1)*3+j);
        %dscatter(dgMat(j,GP)',fccMat((j-1)*3+i,GP)'); %FCC(J_i,v_j)
        dscatter(dgMax(GP)',fccMat((j-1)*3+i,GP)');
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
        dscatter(dgMat(j,GP)',fccMat((j-1)*3+i,GP)'); %FCC(J_i,v_j)
        ylim([-2 2]);
        title(strcat('C_{v_',num2str(j),'}^{J_',num2str(i),'}'));
        xlabel(strcat('g_',num2str(j)));
        ylabel('FCC');
        box on;
    end
end
colormap('redbluecmap');

figure;
for i=1:3
    for j=1:3
        subplot(3,3,(i-1)*3+j);
        SubFCC=fccMat((j-1)*3+i,GP);
        xmin=quantile(SubFCC,0.05);
        xmax=quantile(SubFCC,0.95);
        if xmax<0
            xmax=0;
        elseif xmin>0
            xmin=0;
        end
        if xmin==0 && xmax<1
            xmax=1;
        end
        hist(SubFCC(SubFCC>=xmin & SubFCC<=xmax),xmin+(xmax-xmin)/40:...
            (xmax-xmin)/20:xmax-(xmax-xmin)/40); 
        xlim([xmin xmax]);
        title(strcat('FCC(',num2str(i),',',num2str(j),')'));
    end
end