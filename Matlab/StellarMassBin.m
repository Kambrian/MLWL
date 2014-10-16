figure;
[xm1,ym1]=linhist(gal.Zspec(gal.SM<=10.5),50,'stairs','r');
hold on;
[xm2,ym2]=linhist(gal.Zspec(gal.SM>10.5&gal.SM~=9999),50,'stairs','g');
dlmwrite('GalBinZdistr.dat',[xm1,ym1,xm2,ym2]);

%%
myfigure;
[h,r1,s1,es1,cov1]=plotgama_prof('GALlow.ZGAMA.V4','r.');
hold on;
[h,r2,s2,es2,cov2]=plotgama_prof('GALhigh.ZGAMA.V4','go');
%%
dir='/home/kam/Projects/Lensing/data/stellarmass_data/fig1/';
% file=[dir,'rebin.sm.all.sm4.highfdev.dat'];
% data=importdata(file);

% figure;loglog(data(:,1)/1000,data(:,2)*100);
% hold on;
file=[dir,'fitavgsig.hh.all.sm4.highfdev.dat'];
data=importdata(file);
data=[data(2:end,:);data(1,:)];
semilogx(data(:,1),log10(data(:,2)*100),'r-');
hold on;

file=[dir,'fitavgsig.hh.all.sm4.lowfdev.dat'];
data=importdata(file);
data=[data(2:end,:);data(1,:)];
loglog(data(:,1),log10(data(:,2)*100),'r--');
hold on;

file=[dir,'fitavgsig.hh.all.sm5.highfdev.dat'];
data=importdata(file);
data=[data(2:end,:);data(1,:)];
loglog(data(:,1),log10(data(:,2)*100),'g-');
hold on;

file=[dir,'fitavgsig.hh.all.sm5.lowfdev.dat'];
data=importdata(file);
data=[data(2:end,:);data(1,:)];
loglog(data(:,1),log10(data(:,2)*100),'g--');
hold on;

file=[dir,'fitavgsig.hh.all.sm6.highfdev.dat'];
data=importdata(file);
data=[data(2:end,:);data(1,:)];
loglog(data(:,1),log10(data(:,2)*100),'b-');
hold on;

file=[dir,'fitavgsig.hh.all.sm6.lowfdev.dat'];
data=importdata(file);
data=[data(2:end,:);data(1,:)];
loglog(data(:,1),log10(data(:,2)*100),'b--');
hold on;

% legend('M*<10.5','M*>10.5')
print('-depsc','GalLensSMbin.eps');