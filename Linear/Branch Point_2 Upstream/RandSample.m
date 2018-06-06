NSample=20000;

mu=1;
sigma=1;
k=10.^normrnd(mu,sigma,3,NSample);
K=10.^normrnd(mu,sigma,3,NSample);
Sin1=10.^normrnd(mu,sigma,1,NSample);
Sin2=10.^normrnd(mu,sigma,1,NSample);
Sout=10.^normrnd(mu,sigma,1,NSample);

for i=1:NSample
    [fcc,J,dg]=MCA_Converged(k(:,i),K(:,i),Sin1(i),Sin2(i),Sout(i));
    fccMat(:,i)=fcc(:);
    JMat(:,i)=J(:);
    dgMat(:,i)=dg(:);
end

dgMax=max(dgMat);
dgMin=min(dgMat);
GP=find(min(JMat)>0); %GoodPos, which means fluxes through each reaction is positive
figure;
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
        if abs(xmax-1)<0.09
            xmax=1;
        end
        hist(SubFCC(SubFCC>=xmin & SubFCC<=xmax),xmin+(xmax-xmin)/40:...
            (xmax-xmin)/20:xmax-(xmax-xmin)/40); 
        xlim([xmin xmax]);
        title(strcat('FCC(',num2str(i),',',num2str(j),')'));
    end
end

figure;
for i=1:3
    for j=1:3
        subplot(3,3,(i-1)*3+j);
        dscatter(log10(JMat(j,GP)'),fccMat((j-1)*3+i,GP)'); %FCC(J_i,v_j)
        ylim([-1 1]);
        title(strcat('FCC(',num2str(i),',',num2str(j),')'));
        xlabel(strcat('J_',num2str(j)));
        ylabel('FCC');
    end
end
colormap('redbluecmap');

