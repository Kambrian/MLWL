load('additional/database/smhm_guo13.mat');
load('additional/database/smhm_bower.mat');
load('additional/database/smhm_font.milli.mat');
%%
file='additional/galform/lagos/gal_vol3.hdf5';
mstar=h5read(file,'/Output001/mstars_tot');
mhalo=h5read(file,'/Output001/mhhalo');
SFR=h5read(file,'/Output001/mstardot');
iscen=h5read(file,'/Output001/is_central');
f=iscen>0;
lagos.Mstar=mstar(f)*0.73; %convert to Msun/h^2
lagos.Mhalo=mhalo(f); %Msun/h
lagos.SFR=SFR(f)*0.73; %Msun/h^2/Gyr
lagos.SSFR=lagos.SFR./lagos.Mstar; %1/Gyr
%%
clear
clc
addpath(genpath('/work/Projects/Lensing/code/v8.0/Matlab'))
cd /work/Projects/Lensing/data
load G3Cv4up8/mockcat_1.mat
load G3Cv4up8/G3Cv4up8.mat
%%
prepare_smhm_HODs;
%% Success
data=[7.3685e+09   1.45222e+09  11.51      0.5988                  
1.4683e+10   2.90974e+09  11.97      0.3629                  
2.57841e+10   3.59833e+09  12.16      0.3087                  
4.16434e+10   5.83513e+09  12.43      0.2122                  
6.67911e+10   9.34632e+09  12.61      0.2369                  
1.07375e+11   1.52612e+10  13.26      0.1444                  
1.73309e+11   2.39779e+10  13.34      0.2277                  
2.77191e+11   3.80254e+10  nan      nan                   
4.48464e+11   5.69671e+10  12.9       1.273
];%r<19.4 centrals
xbrd=[5e9,1e10; 1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11; 3.7e11, 6e11];
xbrd=log10(xbrd);
datahigh=[7.29997e+09   1.43735e+09  12.01      0.4588                  
1.442e+10   2.88151e+09  12.51      0.3106                  
2.52676e+10   3.51158e+09  12.8       0.3169                  
3.97499e+10   5.5919e+09  13.11      0.3219                  
6.31942e+10   8.748e+09  5.297     nan                   %266 halos
1.00395e+11   1.35373e+10  13         nan                   %41 halos
1.5887e+11   1.47669e+10  14.22      1.096                   %3 halos in total
2.61078e+11   2.69671e+10  10.28      nan                    %3 halos in total
nan nan nan nan
]; %high SSFR, and restricted to r<19.4
datacomp=[7.44261e+09   1.46453e+09  8.88       nan                   
1.48636e+10   2.91526e+09  11.35      1.267                   
2.59596e+10   3.61056e+09  11.96      0.4464                  
4.18871e+10   5.82114e+09  12.35      0.2438                  
6.6924e+10   9.3412e+09  12.64      0.23                    
1.07442e+11   1.52613e+10  13.25      0.1453                  
1.73332e+11   2.3983e+10  13.34      0.2288                  
2.77271e+11   3.80553e+10  5.082     nan                   
4.48464e+11   5.69671e+10  12.9       1.273 ]; %complementary to high SSFR, within r<19.4 centrals
myfigure;
sigma_common=0.2;
for i=1:numel(ms)
    zwang(i)=medianMh(log10(ms(i)),sigma_common,Pwang4);
    zling(i)=medianMh(log10(ms(i)),sigma_common,Pling);
    zmosterZ0(i)=medianMh(log10(ms(i)),sigma_common,PmosterZ0);
    zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo); 
end
zbound=minmax([zwang',zling',zmosterZ0',zguo']);
h3=area(log10(ms),[zbound(:,1),zbound(:,2)-zbound(:,1)]);hold on;
set(h3(2),'facecolor','y','edgecolor','w');
set(h3(1),'facecolor','w','edgecolor','w');
set(gca,'layer','top');
%----------mock predictions----------------
x=linspace(8,12,12);
%-------standard mock--------------
f=galmock.CentralSampleIter>0; %about 90% are real centrals
[xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(galmock.Mstar(f)),log10(galmock.Mhalo(f)),x,0.68);
h2=plot(xmed,ymed,'-','linewidth',2,'color','r');hold on;
plot(xmed,yl,'r--');

h1=ploterr(log10(data(:,1)),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro');hold on;
% f=datahigh(:,1)<1e12;
% h11=ploterr(log10(datahigh(f,1)),datahigh(f,3),{xbrd(f,1),xbrd(f,2)},datahigh(f,4),'bs','logx');hold on;
xlim([9.5,12]);
ylim([10,15.5]);
% xscale('log');
xlabel('$\log$(Stellar Mass[M$_\odot/h^2$])','fontsize',25);
ylabel('$\log$(Halo Mass[M$_\odot/h$])','fontsize',25);
l=legend([h1(1),h2, h3(2)],'GAMA Central','Mock Central','$\sigma=0.2$ HOD');
set(l,'location','northwest','interpreter','latex','fontsize',25);
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/presentation/M-Mstar-central.eps');
%% HOD
%% real central
data=[7.3685e+09   1.45222e+09  11.51      0.5988                  
1.4683e+10   2.90974e+09  11.97      0.3629                  
2.57841e+10   3.59833e+09  12.16      0.3087                  
4.16434e+10   5.83513e+09  12.43      0.2122                  
6.67911e+10   9.34632e+09  12.61      0.2369                  
1.07375e+11   1.52612e+10  13.26      0.1444                  
1.73309e+11   2.39779e+10  13.34      0.2277                  
2.77191e+11   3.80254e+10  nan      nan                   
4.48464e+11   5.69671e+10  12.9       1.273
];%r<19.4 centrals
xbrd=[5e9,1e10; 1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11; 3.7e11, 6e11];
xbrd=log10(xbrd);
datahigh=[7.29997e+09   1.43735e+09  12.01      0.4588                  
1.442e+10   2.88151e+09  12.51      0.3106                  
2.52676e+10   3.51158e+09  12.8       0.3169                  
3.97499e+10   5.5919e+09  13.11      0.3219                  
6.31942e+10   8.748e+09  5.297     nan                   %266 halos
1.00395e+11   1.35373e+10  13         nan                   %41 halos
1.5887e+11   1.47669e+10  14.22      1.096                   %3 halos in total
2.61078e+11   2.69671e+10  10.28      nan                    %3 halos in total
nan nan nan nan
]; %high SSFR, and restricted to r<19.4
datacomp=[7.44261e+09   1.46453e+09  8.88       nan                   
1.48636e+10   2.91526e+09  11.35      1.267                   
2.59596e+10   3.61056e+09  11.96      0.4464                  
4.18871e+10   5.82114e+09  12.35      0.2438                  
6.6924e+10   9.3412e+09  12.64      0.23                    
1.07442e+11   1.52613e+10  13.25      0.1453                  
1.73332e+11   2.3983e+10  13.34      0.2288                  
2.77271e+11   3.80553e+10  5.082     nan                   
4.48464e+11   5.69671e+10  12.9       1.273 ]; %complementary to high SSFR, within r<19.4 centrals
datahighsat=[7.1879e+09   1.45482e+09  6.042      nan%308                     
1.41622e+10   2.81735e+09  12.48      0.5461                  
2.52305e+10   3.54182e+09  12.44      0.79                    
3.873e+10   4.94786e+09  12.82      1.013                   
6.06589e+10   6.96067e+09  12.97      1.647                   
9.94397e+10   1.11144e+10  13.89      1.284                   
1.58129e+11   0  15.48      0.8361
nan nan nan nan
nan nan nan nan]; %high SSFR satellites, r<19.4
datasat=[7.33134e+09   1.43807e+09  5.807      nan%197                     
1.46576e+10   2.85986e+09  12.7       0.1572                  
2.5628e+10   3.56798e+09  12.55      0.2525                  
4.15879e+10   5.79792e+09  13.12      0.1474                  
6.60348e+10   9.32018e+09  12.86      0.27                    
1.05665e+11   1.4216e+10  13.55      0.1878                  
1.68307e+11   2.34732e+10  13.74      0.3021                  
2.71853e+11   3.56047e+10  8.206      nan%51.84                   
4.53185e+11   6.57651e+10  13.73      1.337]; %satellites, r<19.4

myfigure;
% sigma_common=0.2;
% for i=1:numel(ms)
%     zwang(i)=medianMh(log10(ms(i)),sigma_common,Pwang4);
%     zling(i)=medianMh(log10(ms(i)),sigma_common,Pling);
%     zmosterZ0(i)=medianMh(log10(ms(i)),sigma_common,PmosterZ0);
%     zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo); 
% end
% zbound=minmax([zwang',zling',zmosterZ0',zguo']);
% h3=area(log10(ms),[zbound(:,1),zbound(:,2)-zbound(:,1)]);hold on;
% set(h3(2),'facecolor','y','edgecolor','w');
% set(h3(1),'facecolor','w','edgecolor','w');
% set(gca,'layer','top');
%----------mock predictions----------------
x=linspace(8,12,12);
%---real centrals----- %little difference from the CentralSampleIter
% f=galmock.is_central>0;
% [xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(galmock.Mstar(f)),log10(galmock.Mhalo(f)),x,0.68);
% h22=plot(xmed,ymed,'--','linewidth',2,'color','r');hold on;
%-------standard mock--------------
f=galmock.CentralSampleIter>=0&galmock.is_central==1; %about 90% are real centrals
[xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(galmock.Mstar(f)),log10(galmock.Mhalo(f)),x,0.68);
h3=plot(xmed,ymed,'--','linewidth',4,'color','r');hold on;
% plot(xmed,yl,'r--');
f=galmock.CentralSampleIter>=0&galmock.is_central==0;
[xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(galmock.Mstar(f)),log10(galmock.Mhalo(f)),x,0.68);
h33=plot(xmed,ymed,'--','linewidth',4,'color','b');hold on;
%----------sat
f=galmock.CentralSampleIter>0; %about 90% are real centrals
[xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(galmock.Mstar(f)),log10(galmock.Mhalo(f)),x,0.68);
h2=plot(xmed,ymed,'-','linewidth',2,'color','r');hold on;
% plot(xmed,yl,'r--');
f=galmock.CentralSampleIter==0;
[xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(galmock.Mstar(f)),log10(galmock.Mhalo(f)),x,0.68);
h22=plot(xmed,ymed,'-','linewidth',2,'color','b');hold on;
% plot(xmed,yl,'b--');
% f=(galmock.EnvLevel==0|galmock.IsIterCen>0)&galmock.SFR./galmock.Mstar<10^-1.5;
% x=logspace(8,12,15);
% [xmed,ymed,yl,xm,ym,yerr]=skeleton(galmock.Mstar(f),log10(galmock.Mhalo(f)),x,0.68);
% h222=plot(xmed,ymed,'-','linewidth',2,'color','g');hold on;
% ------------volume limited mock----------------
% tmpmock=bower;
% tmpmock.Mstar=tmpmock.sm;tmpmock.Mhalo=tmpmock.mh;tmpmock.SSFR=tmpmock.sfr./tmpmock.sm;
% [xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(tmpmock.Mstar),log10(tmpmock.Mhalo),x,0.68);
% h2=plot(xmed,ymed,'--','linewidth',4,'color','r');hold on;
% % plot(xmed,yl,'r--');
% f=tmpmock.SSFR>10^-1.5;
% [xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(tmpmock.Mstar(f)),log10(tmpmock.Mhalo(f)),x,0.68);
% h22=plot(xmed,ymed,'--','linewidth',4,'color','b');hold on;
% plot(xmed,yl,'b--');
%---------HOD models---------------------
% plot(xmed,13.55+0.81*log10(xmed/1e12),'k-');
% plot(xmed,13.85+0.92*log10(xmed/1e12),'k--');
% h4=plot(ms,ywang4,'c--');
% h5=plot(ms,yling,'m--');
% h6=plot(ms,ymosterZ0,'g--');
% mh=logspace(10,16,10);
% h7=plot(log10(halo2starP(mh,Pguo)),log10(mh),'b--'); 


h1=ploterr(log10(data(:,1)),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro');hold on;
h11=ploterr(log10(datasat(:,1)),datasat(:,3),{xbrd(:,1),xbrd(:,2)},datasat(:,4),'co');hold on;
% h11=ploterr(log10(datahighsat(:,1)),datahighsat(:,3),{xbrd(:,1),xbrd(:,2)},datahighsat(:,4),'gd','logx');hold on;
f=datahigh(:,1)<1e12;
h111=ploterr(log10(datahigh(f,1)),datahigh(f,3),{xbrd(f,1),xbrd(f,2)},datahigh(f,4),'bs','logx');hold on;
% h111=ploterr(log10(dataactive(f,1)),dataactive(f,3),{xbrd(f,1),xbrd(f,2)},dataactive(f,4),'mo','logx');hold on;
% h=ploterr(log10(datahighz(:,1)),datahighz(:,3),[],[],'bx');set(h(1),'markersize',15);
% plot(x,15.2+1.58*(x-12),'b-');
% plot(x,15.01+1.43*(x-12),'r-');
% plot(x,13.96+0.93*(x-12),'g-');
% plot(x,13.75+0.78*(x-12),'g--');
% plot(x,13.59+0.85*(x-12),'r-');
xlim([9.5,12]);
ylim([10,15.5]);
% xscale('log');
xlabel('$\log(M_\star$[M$_\odot/h^2$])');
ylabel('$\log(M_h$[M$_\odot/h$])');
l=legend([h1(1),h11(1),h111(1),h2,h22,h3,h33],'All Central','All Sat','Active Central','Mock central','Mock Sat','Mock true central','Mock true sat');
% l=legend([h1(1),h2,h22,h222,h4,h5,h6,h7,h3(2)],'Data','Mock','Mock Active','Mock Passive','WangL13','WangLY13','Moster13','Guo10','All $\sigma=0.2$');
set(l,'location','northwest','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/presentation/M-Mstar-sat.eps');
%% scatter/contour plot
nbin=50;
smth=@(n) conv2(n,[0.05,0.1,0.05;0.1,0.4,0.1;0.05,0.1,0.05;],'same');
percents=[0.3,0.6,0.9];
myfigure;
f=galmock.CentralSampleIter>=0&galmock.SSFR<10^-1.5&galmock.is_central==0;
x=log10(galmock.Mstar(f));y=log10(galmock.Mhalo(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[9.5,12],[10,15.5]);
[~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'m','linewidth',1);hold on;%colormap(cmap);
f=galmock.CentralSampleIter>=0&galmock.SSFR>10^-1.5&galmock.is_central==0;
x=log10(galmock.Mstar(f));y=log10(galmock.Mhalo(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[9.5,12],[10,15.5]);
[~,h2]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'g','linewidth',1);hold on;%colormap(cmap);
x=linspace(8,12,12);
f=galmock.CentralSampleIter>=0; %about 90% are real centrals
[xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(galmock.Mstar(f)),log10(galmock.Mhalo(f)),x,0.68);
h3=plot(xmed,ymed,'-','linewidth',2,'color','r');hold on;
h11=ploterr(log10(data(:,1)),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro');hold on;
f=datahigh(:,1)<1e12;
h111=ploterr(log10(datahigh(f,1)),datahigh(f,3),{xbrd(f,1),xbrd(f,2)},datahigh(f,4),'bs','logx');hold on;
l=legend([h111(1),h11(1),h2,h1,h3],'Active lensing','All lensing','SatActive mock','SatPassive mock','All mock(median)');set(l,'fontsize',20,'location','south','box','off');
xlim([9.5,12]);ylim([10,15.5]);
xlabel('$\log($Stellar Mass[M$_\odot/h^2$])','fontsize',25);
ylabel('$\log$(Halo Mass[M$_\odot/h$])','fontsize',25);
% print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/presentation/M-Mstar-contour-sat.eps');
%%
nbin=30;
smth=@(n) conv2(n,[0.05,0.1,0.05;0.1,0.4,0.1;0.05,0.1,0.05;],'same');
% smth=@(n) conv2(n,[0.08,0.15,0.08;0.15,0.23,0.15;0.08,0.15,0.08;],'same');
percents=[0.3,0.6,0.9];
myfigure;
f=gal.FlagSFR>0;
x=log10(gal.SMsps(f));y=log10(gal.SFR(f)./gal.SMsps(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[-5,4]);
[~,h0]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'r','linewidth',2);hold on;%colormap(cmap);
f=galmock.zhub<0.31;
x=log10(galmock.Mstar(f));y=log10(galmock.SFR(f)./galmock.Mstar(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[-5,4]);
[~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'g','linewidth',2);hold on;%colormap(cmap);
% f=rawmock.z<0.31;
% x=log10(rawmock.sm(f));y=log10(rawmock.SFR(f)./rawmock.sm(f));
% [xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[-5,4]);
% [~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'b');hold on;%colormap(cmap);
plot([8,12],[-1.5,-1.5],'k--');
l=legend('Data','Mock');set(l,'interpreter','latex','fontsize',25);
xlabel('$\log(M_\star$[M$_\odot/h^2$])','fontsize',25);
ylabel('$\log$(SSFR[Gyr$^{\frac{\ }{\ }1}$])','fontsize',25);
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/presentation/SSFR-Mstar.eps');