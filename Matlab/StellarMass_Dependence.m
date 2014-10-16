clear
clc
addpath(genpath('/work/Projects/Lensing/code/v8.0/Matlab'))
cd /work/Projects/Lensing/data
stars=[];
for vol=1:9
    load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
    f=grpmock.Mult>2;
    stars=[stars;grpmock.MstarIter(f),grpmock.MIter(f),grpmock.SFRIter(f)];
end
ssfr=stars(:,3)./stars(:,1);
load G3Cv4up8/mockcat_1.mat
load G3Cv4up8/G3Cv4up8.mat
%%
bower=importdata('smhm_bower.dat');
bower.LumHalpha=bower.data(:,4)*1e33; %J/s/h2
bower.sfr=bower.LumHalpha/1.27e34*1e9; %Msun/h2/Gyr
bower.sm=bower.data(:,3)*0.7; %Msun/h2
bower.mh=bower.data(:,1)*8.6e8; %Msun/h
bower.md=bower.data(:,2); %Msun/h, the dhalo mass used in the galform model, required to be non-decreasing so is no smaller than mh.
%%
data=importdata('/work/Projects/Lensing/code/v6.0/Matlab/cosmology/dndlogM_ST_3778.dat',' ',1);
hmf=data.data;
sp=spline(log10(hmf(:,1)),log(hmf(:,3)));
halomfunc=@(logm) exp(ppval(sp,logm)); %dn/dlog10(M)/dV, for halo mass func
clear data
newMF=load('/work/Projects/Lensing/data/StellarMassFunc/massfun.dr72bbright0.modelmass');
%
MvInMvb=0.893;
MvcInMvb=0.733;
halo2star=@(M,A,M0,alpha,beta,gamma) A./((M/M0).^alpha+(M/M0).^beta).^gamma.*M;
halo2starP=@(M,P) halo2star(M,P(1),P(2),P(3),P(4),P(5)); %Msun/h^2 for Mstar, Msun/h for Mhalo
dlnMhdlnMs=@(Mh,P) 1./(1-P(5).*(P(3)*(Mh/P(2)).^P(3)+P(4)*(Mh/P(2)).^P(4))./((Mh/P(2)).^P(3)+(Mh/P(2)).^P(4)));
Pwang4=[2*10^10.23/3.43e11*0.73^2*MvcInMvb,3.43e11/MvcInMvb,1-2.56,1-0.34,1];%DR7, 2013
Pguo=[0.129*0.73*MvcInMvb,10^11.4*0.73/MvcInMvb,-0.926,0.261,2.440]; %close to Yang at >M0 ; M200c
Pling=[2*10^-1.73,10^11.70,-1.16,0.71,1];%M200b
z=0.;
parmosterZ=[0.0351-0.0247*z/(1+z),11.590+1.195*z/(1+z),1.376-0.826*z/(1+z),0.608+0.329*z/(1+z)];
PmosterZ0=[2*parmosterZ(1)*0.72*MvcInMvb,10^parmosterZ(2)*0.72/MvcInMvb,-parmosterZ(3),parmosterZ(4),1];%Moster13, M200c
%
powerlaw_err=@(x,lnAerr,Berr,corrAB) sqrt(lnAerr^2+(log(x)*Berr).^2+2*log(x)*corrAB*lnAerr*Berr);
powerlawmodel_err=@(m,perr,c) powerlaw_err(m/1e12,perr(1)/log10(exp(1)),perr(2),c); %error in ln(y) for y=10^p(1)*(m/1e12).^p(2)
powerlawmodel=@(m,p) 10^p(1)*(m/1e12).^p(2);
%
logMmin=7;
logMmax=17;
sigma0=0.0;%additional scatter due to our SM estimation method
dndlogMs=@(logms,sigma,par) quadgk(@(logm) halomfunc(logm).*normpdf(logms-log10(halo2starP(10.^logm,par)),0,sqrt(sigma^2+sigma0^2)),logMmin,logMmax);
dndlogMsdlogMh=@(logms,logmh,sigma,par) halomfunc(logmh).*normpdf(logms-log10(halo2starP(10.^logmh,par)),0,sqrt(sigma^2+sigma0^2));
medianMh=@(logMs,sigma,par) quadgk(@(logmh) dndlogMsdlogMh(logMs,logmh,sigma,par)./dndlogMs(logMs,sigma,par).*logmh, logMmin, logMmax);
halo2smfunc=@(mh,par) halomfunc(log10(mh)).*abs(dlnMhdlnMs(mh,par));

ms=logspace(8,12,20);
ywang4=ms;sigma_wang4=0.17; %Wang13
yling=ms;sigma_ling=0.22;
ymosterZ0=ms;sigma_mosterZ0=0.1;
yguo=ms;sigma_guo=0.0; 
for i=1:numel(ms)
    ywang4(i)=medianMh(log10(ms(i)),sigma_wang4,Pwang4);
    yling(i)=medianMh(log10(ms(i)),sigma_ling,Pling);
    ymosterZ0(i)=medianMh(log10(ms(i)),sigma_mosterZ0,PmosterZ0);
    yguo(i)=medianMh(log10(ms(i)),sigma_guo,Pguo); 
end
%% common dispersion
sigma_common=0.1;
zwang=ms;
zling=ms;
zmosterZ0=ms;
zguo=ms;
for i=1:numel(ms)
    zwang(i)=medianMh(log10(ms(i)),sigma_common,Pwang4);
    zling(i)=medianMh(log10(ms(i)),sigma_common,Pling);
    zmosterZ0(i)=medianMh(log10(ms(i)),sigma_common,PmosterZ0);
    zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo); 
end
zbound=minmax([zwang',zling',zmosterZ0',zguo']);
%% M=A(M*)
data=[1.51429e+10   2.91794e+09  12.24      1.216
    2.61354e+10   3.59945e+09  12.15      1.12                    
4.24234e+10   5.81214e+09  12.76      0.3252                  
6.77309e+10   9.51513e+09  13.1       0.227                   
1.08512e+11   1.54176e+10  13.63      0.1394                  
1.73428e+11   2.33892e+10  13.58      0.2571                       
2.7515e+11   3.71434e+10  14.09      0.2438                  
% 4.53403e+11   5.96044e+10  9.578      51.67 
]; %sps mass
% data=[2.62242e+10   3.49365e+09  12.45      0.5894                  
% 4.20769e+10   5.80834e+09  12.76      0.3103                  
% 6.73279e+10   9.4463e+09  13.32      0.1708                  
% 1.09286e+11   1.48898e+10  13.65      0.1524                  
% 1.71909e+11   2.44146e+10  13.75      0.2274                  
% 2.79703e+11   3.79262e+10  13.59      0.6948                  
% 4.39652e+11   4.87496e+10  14.41      0.6111]; %color mass 
% data=[2.60354e+10   3.60113e+09  12.39      0.4289                  
% 4.17511e+10   5.84704e+09  12.59      0.288                   
% 6.70784e+10   9.42636e+09  13.04      0.1992                  
% 1.0811e+11   1.48238e+10  13.48      0.1614                  
% 1.71039e+11   2.37215e+10  13.25      0.4205                  
% 2.7939e+11   3.79044e+10  13.37      0.807                   
% 4.4281e+11   5.69634e+10  13.99      0.8274]; %N>1
data2=[%  0 |     a =  11.74   |  0.2142  |          |          |          |          |          |
    1.51429e+10   2.91794e+09  11.67      2.096 %first bin is all the central galaxies (down to N=1) not included in other bins.
    2.61354e+10   3.59945e+09  7.206      nan                    
4.24234e+10   5.81214e+09  12.87    0.2828                 
6.77309e+10   9.51513e+09  12.98    0.272                   
1.08512e+11   1.54176e+10  13.65    0.1387                  
1.73428e+11   2.33892e+10  13.49   0.274                       
2.7515e+11   3.71434e+10  13.77      0.4103 ]; %jointly fitted bins
xbrd=[1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11];
myfigure;
for vol=1 %:9
    load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
f=grpmock.Mult>2;
x=logspace(9.8,12,10);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.MstarIter(f),log10(grpmock.MIter(f)),x,0.68);
h3=plot(xmed,ymed,'g');hold on;
end
plot(xmed,yl,'g--');
h1=ploterr(data(:,1),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro','logx');hold on;
% h2=plot(xmed,15.06+1.57*log10(xmed/1e12),'k-'); %color mass
% h2=plot(xmed,14.42+1.08*log10(xmed/1e12),'k-');%sps mass
% h22=plot(xmed,14.42-0.05+1.1*log10(xmed/1e12),'k--');%sps mass
h11=ploterr(data2(:,1),data2(:,3),[],data2(:,4),'bs','logx');hold on;
xlim([1e10,1e12]);
ylim([10.8,15.5]);
xscale('log');
xlabel('$M_\star$[M$_\odot/h^2$]');
ylabel('$\log(M_h$[M$_\odot/h$])');
l=legend([h1(1),h11(1)],'Binned','Joint Binned','Systematic','Mock');set(l,'location','southeast','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/extra/M-Mstar-jointbin.eps');
%% M=A(M*) comparison
data=[1.51429e+10   2.91794e+09  12.24      1.216
    2.61354e+10   3.59945e+09  12.15      1.12                    
4.24234e+10   5.81214e+09  12.76      0.3252                  
6.77309e+10   9.51513e+09  13.1       0.227                   
1.08512e+11   1.54176e+10  13.63      0.1394                  
1.73428e+11   2.33892e+10  13.58      0.2571                       
2.7515e+11   3.71434e+10  14.09      0.2438                  
% 4.53403e+11   5.96044e+10  9.578      51.67 
]; %sps mass
xbrd=[1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11];
myfigure;
f=grpmock.Mult>2;
x=logspace(9,12,15);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.MstarIter(f),log10(grpmock.MIter(f)),x,0.68);
h2=plot(xmed,ymed,'-','linewidth',2,'color','r');hold on;

h4=plot(ms,ywang4,'c--');
h5=plot(ms,yling,'--','color',[0.8,0.8,1]);
h6=plot(ms,ymosterZ0,'g--');
mh=logspace(10,16,10);
h7=plot(halo2starP(mh,Pguo),log10(mh),'b--'); 
% sigma_common=0.1;
% for i=1:numel(ms)
%     zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo);
% end
% plot(ms,zguo,'g:');
% sigma_common=0.2;
% for i=1:numel(ms)
%     zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo);
% end
% plot(ms,zguo,'g-.');
% sigma_common=0.3;
% for i=1:numel(ms)
%     zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo);
% end
% plot(ms,zguo,'g-');

h1=ploterr(data(:,1),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro','logx');hold on;
% set(h1(1),'markersize',7)
xlim([9e9,1e12]);
ylim([10.8,15.5]);
xscale('log');
xlabel('$M_\star$[M$_\odot/h^2$]');
ylabel('$\log(M_h$[M$_\odot/h$])');
% xlabel('Central Stellar Mass[M$_\odot/h^2$]');
% ylabel('$\log($Halo Mass[M$_\odot/h$])');
l=legend([h1(1),h2,h7,h6,h4,h5],'Data','Mock','Guo10','Moster13','WangL13','WangLY13');
set(l,'location','southeast','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/M-Mstar-comp.eps');
%% M=A(M*) comparison (common dispersion)
data=[1.51429e+10   2.91794e+09  12.24      1.216
    2.61354e+10   3.59945e+09  12.15      1.12                    
4.24234e+10   5.81214e+09  12.76      0.3252                  
6.77309e+10   9.51513e+09  13.1       0.227                   
1.08512e+11   1.54176e+10  13.63      0.1394                  
1.73428e+11   2.33892e+10  13.58      0.2571                       
2.7515e+11   3.71434e+10  14.09      0.2438                  
];
data5=[1.52189e+10   2.77664e+09  12.92      1.244                   
2.64628e+10   3.86419e+09  12.2       1.607                   
4.41984e+10   5.89932e+09  13.11      0.4577                  
6.90324e+10   9.4982e+09  13.44      0.2624                  
1.09594e+11   1.60071e+10  13.88      0.1412                  
1.74127e+11   2.34103e+10  13.99      0.1942                  
2.73782e+11   3.49934e+10  14.36      0.2142                  
% 4.33703e+11   4.69153e+10  9.908      63.45 
]; %N>5;
xbrd=[1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11];
myfigure;
% sigma_common=0.1;
% for i=1:numel(ms)
%     zwang(i)=medianMh(log10(ms(i)),sigma_common,Pwang4);
%     zling(i)=medianMh(log10(ms(i)),sigma_common,Pling);
%     zmosterZ0(i)=medianMh(log10(ms(i)),sigma_common,PmosterZ0);
%     zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo); 
% end
% zbound=minmax([zwang',zling',zmosterZ0',zguo']);
% h8=area(ms,[zbound(:,1),zbound(:,2)-zbound(:,1)]);hold on;
% set(h8(2),'facecolor','c','edgecolor','w');
% set(h8(1),'facecolor','w','edgecolor','w');
sigma_common=0.2;
for i=1:numel(ms)
    zwang(i)=medianMh(log10(ms(i)),sigma_common,Pwang4);
    zling(i)=medianMh(log10(ms(i)),sigma_common,Pling);
    zmosterZ0(i)=medianMh(log10(ms(i)),sigma_common,PmosterZ0);
    zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo); 
end
zbound=minmax([zwang',zling',zmosterZ0',zguo']);
h3=area(ms,[zbound(:,1),zbound(:,2)-zbound(:,1)]);hold on;
set(h3(2),'facecolor','y','edgecolor','w');
set(h3(1),'facecolor','w','edgecolor','w');
sigma_common=0.35;
for i=1:numel(ms)
    zwang(i)=medianMh(log10(ms(i)),sigma_common,Pwang4);
    zling(i)=medianMh(log10(ms(i)),sigma_common,Pling);
    zmosterZ0(i)=medianMh(log10(ms(i)),sigma_common,PmosterZ0);
    zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo); 
end
zbound=minmax([zwang',zling',zmosterZ0',zguo']);
h4=area(ms,[zbound(:,1),zbound(:,2)-zbound(:,1)]);hold on;
set(h4(2),'facecolor','b','edgecolor','w');
set(h4(1),'facecolor','w','edgecolor','w');
set(gca,'layer','top');
f=grpmock.Mult>5;
x=logspace(9.8,12,15);
[xmed,ymed]=skeleton(grpmock.MstarIter(f),log10(grpmock.MIter(f)),x,0.68);
h00=plot(xmed,ymed,'-','linewidth',2,'color','g');hold on;
f=grpmock.Mult>2;
x=logspace(9.8,12,15);
[xmed,ymed]=skeleton(grpmock.MstarIter(f),log10(grpmock.MIter(f)),x,0.68);
h11=plot(xmed,ymed,'-','linewidth',2,'color','r');hold on;
f=galmock.is_central>0;
% plot(galmock.Mstar(f&galmock.SFR./galmock.Mstar<4e-2),log10(galmock.Mhalo(f&galmock.SFR./galmock.Mstar<4e-2)),'g.','markersize',1);
% plot(galmock.Mstar(f&galmock.SFR./galmock.Mstar>4e-2),log10(galmock.Mhalo(f&galmock.SFR./galmock.Mstar>4e-2)),'b.','markersize',1);
% plot(stars(ssfr<1e-2,1),log10(stars(ssfr<1e-2,2)),'b.','markersize',1);
% plot(stars(ssfr>1e-2,1),log10(stars(ssfr>1e-2,2)),'g.','markersize',1);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(galmock.Mstar(f),log10(galmock.Mhalo(f)),x,0.68);
h2=plot(xmed,ymed,'k--');
h0=ploterr(data5(:,1),data5(:,3),[],data5(:,4),'g^','logx');hold on;
h1=ploterr(data(:,1),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro','logx');hold on;
xlim([9e9,1e12]);
ylim([10.8,15.5]);
xscale('log');
xlabel('$M_\star$[M$_\odot/h^2$]');
ylabel('$\log(M_h$[M$_\odot/h$])');
l=legend([h0(1),h1(1),h00,h11,h2,h3(2),h4(2)],'$N>5$','$N>2$','Mock $N>5$','Mock $N>2$','Mock $N>0$','Models $\sigma=0.2$','Models $\sigma=0.35$');
h4=plot(ms,ywang4,'c--');
h5=plot(ms,yling,'m--');
h6=plot(ms,ymosterZ0,'g--');
mh=logspace(10,16,10);
h7=plot(halo2starP(mh,Pguo),log10(mh),'b--'); 

set(l,'location','southeast','interpreter','latex');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/M-Mstar-dispersion.eps');
%% M=A(M*) comparison, GalCat centrals
data0=[1.51429e+10   2.91794e+09  12.24      1.216
    2.61354e+10   3.59945e+09  12.15      1.12                    
4.24234e+10   5.81214e+09  12.76      0.3252                  
6.77309e+10   9.51513e+09  13.1       0.227                   
1.08512e+11   1.54176e+10  13.63      0.1394                  
1.73428e+11   2.33892e+10  13.58      0.2571                       
2.7515e+11   3.71434e+10  14.09      0.2438                  
% 4.53403e+11   5.96044e+10  9.578      51.67 
]; %grpcat 
xbrd0=[1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11];               
data=[7.37402e+09   1.44834e+09  11.54      0.5633                  
1.46631e+10   2.90395e+09  11.9       0.3931                  
2.57748e+10   3.59255e+09  12.18      0.2907                  
4.16293e+10   5.86249e+09  12.47      0.1978                  
6.67569e+10   9.32472e+09  12.65      0.2219                  
1.07459e+11   1.52862e+10  13.23      0.146                   
1.73073e+11   2.38491e+10  13.35      0.221                   
2.76538e+11   3.82676e+10  nan nan                 
4.51297e+11   6.03076e+10  13.26      0.8478
];%galcat, centrals
xbrd=[5e9,1e10; 1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11; 3.7e11, 6e11];
datahigh=[7.29451e+09   1.43144e+09  11.97      0.4627                  
1.43275e+10   2.86376e+09  12.36      0.3665                  
2.51842e+10   3.50482e+09  12.66      0.3725                  
3.96712e+10   5.59058e+09  13.12      0.3165                  
6.32171e+10   8.76871e+09  5.314      nan                   
1.00285e+11   1.33937e+10  13.76      0.7634                  
1.58397e+11   1.28147e+10  nan      nan                   
2.61078e+11   2.69671e+10  nan      nan
nan nan nan nan]; %ssfr>1e-1.5
datacomp=[7.46945e+09   1.46267e+09  10.35      5.001                   
1.48929e+10   2.9089e+09  11.43      1.049                   
2.59527e+10   3.59963e+09  12.05      0.3773                  
4.18519e+10   5.85128e+09  12.4       0.2241                  
6.68703e+10   9.31978e+09  12.68      0.215                   
1.07519e+11   1.52869e+10  13.22      0.1478                  
1.731e+11   2.38564e+10  13.35      0.2103                  
2.76599e+11   3.82935e+10  5.14       nan                   
4.51297e+11   6.03076e+10  13.26      0.8478]; %complementary to ssfr, the remaining part of galcen
datazbound=[7.36217e+09   1.44739e+09  11.53      0.5673                  
1.46087e+10   2.90045e+09  11.85      0.4284                  
2.5687e+10   3.5878e+09  12.21      0.2809                  
4.13266e+10   5.79956e+09  12.38      0.2283                  
6.58923e+10   9.15656e+09  12.63      0.2415                  
1.05298e+11   1.46054e+10  13.36      0.1415                  
1.69995e+11   2.26231e+10  13.54      0.2051                  
2.7353e+11   3.49233e+10  12.65      2.08                    
4.63165e+11   6.3437e+10  13.37      0.9454]; %the same redshift range as the ssfr subsample
myfigure;
sigma_common=0.2;
for i=1:numel(ms)
    zwang(i)=medianMh(log10(ms(i)),sigma_common,Pwang4);
    zling(i)=medianMh(log10(ms(i)),sigma_common,Pling);
    zmosterZ0(i)=medianMh(log10(ms(i)),sigma_common,PmosterZ0);
    zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo); 
end
zbound=minmax([zwang',zling',zmosterZ0',zguo']);
h3=area(ms,[zbound(:,1),zbound(:,2)-zbound(:,1)]);hold on;
set(h3(2),'facecolor','y','edgecolor','w');
set(h3(1),'facecolor','w','edgecolor','w');
set(gca,'layer','top');
%----------mock predictions----------------
% f=galmock.is_central>0;
% x=logspace(8,12,15);
% [xmed,ymed,yl,xm,ym,yerr]=skeleton(galmock.Mstar(f),log10(galmock.Mhalo(f)),x,0.68);
% h22=plot(xmed,ymed,'--','linewidth',2,'color','r');hold on;
for vol=1:9
    load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
f=galmock.EnvLevel==0|galmock.IsIterCen>0; %about 90% are real centrals
x=logspace(8,12,15);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(galmock.Mstar(f),log10(galmock.Mhalo(f)),x,0.68);
h2=plot(xmed,ymed,'-','linewidth',2,'color','r');hold on;
end
% plot(xmed,yl,'r.-');
f=(galmock.EnvLevel==0|galmock.IsIterCen>0)&galmock.SFR./galmock.Mstar>10^-1.5;
x=logspace(8,12,15);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(galmock.Mstar(f),log10(galmock.Mhalo(f)),x,0.68);
h22=plot(xmed,ymed,'-','linewidth',2,'color','b');hold on;
% f=(galmock.EnvLevel==0|galmock.IsIterCen>0)&galmock.SFR./galmock.Mstar<10^-1.5;
% x=logspace(8,12,15);
% [xmed,ymed,yl,xm,ym,yerr]=skeleton(galmock.Mstar(f),log10(galmock.Mhalo(f)),x,0.68);
% h222=plot(xmed,ymed,'-','linewidth',2,'color','g');hold on;
% ---------adding noise in Mstar------------------
galmock.Mstar2=galmock.Mstar.*10.^(normrnd(0,0.3,size(galmock.Mstar)));
galmock.SFR2=galmock.SFR.*10.^(normrnd(0,0.0,size(galmock.Mstar)));
galmock.SFR2=galmock.SFR2(randperm(numel(galmock.Mstar)));
f=galmock.EnvLevel==0|galmock.IsIterCen>0; %about 90% are real centrals
x=logspace(8,12,15);
[xmed,ymed]=skeleton(galmock.Mstar2(f),log10(galmock.Mhalo(f)),x,0.68);
h20=plot(xmed,ymed,'--','linewidth',2,'color','r');hold on;
f=(galmock.EnvLevel==0|galmock.IsIterCen>0)&galmock.SFR2./galmock.Mstar2>10^-1.5;
x=logspace(8,12,15);
[xmed,ymed]=skeleton(galmock.Mstar2(f),log10(galmock.Mhalo(f)),x,0.68);
h220=plot(xmed,ymed,'--','linewidth',2,'color','b');hold on;
% ------------volume limited mock----------------
f=bower.sm>0;
[xmed,ymed]=skeleton(bower.sm(f),log10(bower.mh(f,1)),x,0.68);
h220=plot(xmed,ymed,':','linewidth',2,'color','r');hold on;
f=bower.sfr./bower.sm>10^-1.5;
[xmed,ymed]=skeleton(bower.sm(f),log10(bower.mh(f,1)),x,0.68);
h220=plot(xmed,ymed,':','linewidth',2,'color','b');hold on;
%---------HOD models---------------------
% plot(xmed,13.55+0.81*log10(xmed/1e12),'k-');
% plot(xmed,13.85+0.92*log10(xmed/1e12),'k--');
% h4=plot(ms,ywang4,'c--');
% h5=plot(ms,yling,'m--');
% h6=plot(ms,ymosterZ0,'g--');
% mh=logspace(10,16,10);
% h7=plot(halo2starP(mh,Pguo),log10(mh),'b--'); 


h1=ploterr(data(:,1),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro','logx');hold on;
h11=ploterr(datacomp(:,1),datacomp(:,3),{xbrd(:,1),xbrd(:,2)},datacomp(:,4),'gd','logx');hold on;
% h11=ploterr(datazbound(:,1),datazbound(:,3),{xbrd(:,1),xbrd(:,2)},datazbound(:,4),'gx','logx');hold on;
h111=ploterr(datahigh(:,1),datahigh(:,3),{xbrd(:,1),xbrd(:,2)},datahigh(:,4),'bs','logx');hold on;
xlim([4e9,1e12]);
ylim([10,15.5]);
xscale('log');
xlabel('$M_\star$[M$_\odot/h^2$]');
ylabel('$\log(M_h$[M$_\odot/h$])');
l=legend([h1(1),h111(1),h2,h22,h3(2)],'All Central','Active Central','Mock Central','Mock Active','$\sigma=0.2$ HOD');
% l=legend([h1(1),h2,h22,h222,h4,h5,h6,h7,h3(2)],'Data','Mock','Mock Active','Mock Passive','WangL13','WangLY13','Moster13','Guo10','All $\sigma=0.2$');
set(l,'location','northwest','interpreter','latex');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/M-Mstar-central.eps');
%% compare mock group SSFR dependence
myfigure;
f=grpmock.SSFRIter>10^-1.5; %about 90% are real centrals
x=logspace(8,12,20);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.MstarIter(f),log10(grpmock.Luminosity(f)),x,0.68);
h2=plot(xmed,ymed,'--','linewidth',2,'color','r');hold on;
f=grpmock.Mult>0;
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.MstarIter(f),log10(grpmock.Luminosity(f)),x,0.68);
h2=plot(xmed,ymed,'--','linewidth',2,'color','g');hold on;

f=grp.IterCenSSFR>10^-1.5; %about 90% are real centrals
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grp.IterCenSM(f),log10(grp.Luminosity(f)),x,0.68);
h2=plot(xmed,ymed,'-','linewidth',2,'color','r');hold on;
f=grp.Mult>0;
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grp.IterCenSM(f),log10(grp.Luminosity(f)),x,0.68);
h2=plot(xmed,ymed,'-','linewidth',2,'color','g');hold on;
xscale('log');
%% Nabs dependence in the raw mock, queried from the DGalaxies database, Bower2006a.
data=importdata('Mhalo-Ngal-FixMstar.mock',',',24);
data=data.data;
data(:,1)=data(:,1)*8.6e8;%Msun/h
mbrd=log10(minmax(data(:,3)'));
nbin=logspace(mbrd(1),mbrd(2),40);
colors=spring(numel(nbin));
myfigure;
for i=1:numel(nbin)-1
    f=data(:,3)>=nbin(i)&data(:,3)<nbin(i+1);
loglog(data(f,2),data(f,1),'o','color',colors(i,:));
hold on;
end
xscale('linear');
yscale('log');
xlabel('$M_\star [\rm{M}_\odot/h]$');
ylabel('$M_{h} [\rm{M}_\odot/h]$');
xlim([0.85e11,1.05e11]);
colormap(spring(numel(nbin)));%caxis([1, numel(nbin)]);
hcb = colorbar('YTick',(log10([10,100,1000])-mbrd(1))/(mbrd(2)-mbrd(1))*(numel(nbin)-1)+1,'YTickLabel',...
{'10','100','1000'});
text(10.7e10, 1e11,'$N_{abs}$','fontsize',18)
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/Mhalo-Mstar-N.eps');
%%
myfigure;
% mbrd=(minmax(data(:,2)'));
% nbin=linspace(mbrd(1),mbrd(2),30);
% colors=colormap(jet(numel(nbin)));
% for i=1:numel(nbin)-1
%     f=data(:,2)>=nbin(i)&data(:,2)<nbin(i+1);
% loglog(data(f,3),data(f,1),'o','color',colors(i,:));
% hold on;
% end
loglog(data(:,3),data(:,1),'s');
xscale('log');
yscale('log');
xlabel('$N_{gal}$');
ylabel('$M_{halo} [\rm{M}_\odot/h]$');
% xlim([0.85e11,1.05e11]);
% colormap(jet(numel(nbin)));%caxis([1, numel(nbin)]);
% hcb = colorbar('YTick',(log10([1,10,100,1000])-mbrd(1))/(mbrd(2)-mbrd(1))*(numel(nbin)-1)+1,'YTickLabel',...
% {'N','10','100','1000'});
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3/Mhalo-Ngal-FixMstar.eps');
%%
%% M=A(M*) comparison, for show
data=[1.51429e+10   2.91794e+09  12.24      1.216
    2.61354e+10   3.59945e+09  12.15      1.12                    
4.24234e+10   5.81214e+09  12.76      0.3252                  
6.77309e+10   9.51513e+09  13.1       0.227                   
1.08512e+11   1.54176e+10  13.63      0.1394                  
1.73428e+11   2.33892e+10  13.58      0.2571                       
2.7515e+11   3.71434e+10  14.09      0.2438                  
% 4.53403e+11   5.96044e+10  9.578      51.67 
]; %sps mass
xbrd=[1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11];
myfigure;
f=grpmock.Mult>2;
x=logspace(9.8,12,15);
[xmed,ymed,yl]=skeleton(grpmock.MstarIter(f),grpmock.MIter(f),x,0.68);
h2=loglog(xmed,ymed,'-','linewidth',2,'color','r');hold on;

h4=plot(ms,10.^ywang4,'c--');
h5=plot(ms,10.^yling,'m--');
h6=plot(ms,10.^ymosterZ0,'g--');
mh=logspace(10,16,10);
h7=plot(halo2starP(mh,Pguo),(mh),'b--'); 
% sigma_common=0.1;
% for i=1:numel(ms)
%     zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo);
% end
% plot(ms,zguo,'g:');
% sigma_common=0.2;
% for i=1:numel(ms)
%     zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo);
% end
% plot(ms,zguo,'g-.');
% sigma_common=0.3;
% for i=1:numel(ms)
%     zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo);
% end
% plot(ms,zguo,'g-');
h1=ploterr(data(:,1),10.^data(:,3),{xbrd(:,1),xbrd(:,2)},{10.^(data(:,3)-data(:,4)),10.^(data(:,3)+data(:,4))},'ro','logxy');hold on;
xlim([9e9,1e12]);
ylim(10.^[10.8,15.5]);
xscale('log');yscale('log');
% xlabel('$M_\star$[M$_\odot/h^2$]');
% ylabel('$\log(M_h$[M$_\odot/h$])');
xlabel('Central Stellar Mass[M$_\odot/h^2$]');
ylabel('Halo Mass[M$_\odot/h$])');
l=legend([h1(1),h2,h4,h5,h6,h7],'Data','Mock','WangL13','WangLY13','Moster13','Guo10');
set(l,'location','southeast','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/M-Mstar-comp-show.eps');