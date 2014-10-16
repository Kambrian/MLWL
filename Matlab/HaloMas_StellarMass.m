% dndM=load('/data/jxhan/Downloads/MF_code/test.dndM');
dndM=load('/data/jxhan/Downloads/MF_code/MF3778.dndM');
% dndM=load('/data/jxhan/Downloads/MF_code/millenium.dndM');
sp=spline(log10(dndM(:,1)),log(dndM(:,2).*dndM(:,1).*log(10)));
% data=importdata('/work/Projects/Lensing/code/v6.0/Matlab/cosmology/dndlogM_ST_3778.dat',' ',1);
% data=importdata('/work/Projects/Lensing/code/v6.0/Matlab/cosmology/dndlogM_ST_Planck.dat',' ',1);
% data=importdata('/work/Projects/Lensing/code/v6.0/Matlab/cosmology/dndlogM_ST_WMAP9.dat',' ',1);
% data=importdata('/work/Projects/Lensing/code/v6.0/Matlab/cosmology/dndlogM_ST_WMAP7.dat',' ',1);
% data=importdata('/work/Projects/Lensing/code/v6.0/Matlab/cosmology/dndlogM_ST_Millenium.dat',' ',1);
hmf=data.data;
% sp=spline(log10(hmf(:,1)),log(hmf(:,3)));
halomfunc=@(logm) exp(ppval(sp,logm)); %dn/dlog10(M)/dV, for halo mass func
%%
newMF=load('/work/Projects/Lensing/data/StellarMassFunc/massfun.dr72bbright0.modelmass');
%%
MvInMvb=0.893;
MvcInMvb=0.733;
halo2star=@(M,A,M0,alpha,beta,gamma) A./((M/M0).^alpha+(M/M0).^beta).^gamma.*M;
halo2starP=@(M,P) halo2star(M,P(1),P(2),P(3),P(4),P(5)); %Msun/h^2 for Mstar, Msun/h for Mhalo
dlnMhdlnMs=@(Mh,P) 1./(1-P(5).*(P(3)*(Mh/P(2)).^P(3)+P(4)*(Mh/P(2)).^P(4))./((Mh/P(2)).^P(3)+(Mh/P(2)).^P(4)));
Pwang=[2*10^10.35/4e11*0.73^2*MvcInMvb,4e11/MvcInMvb,1-2.42,1-0.29,1];%2006 quoted in 2010 table ; probably M200c, convert to M200b
Pwang06=[2*10.^10.27/3.33e11*0.73^2*MvcInMvb,3.33e11/MvcInMvb,1-2.59,1-0.276,1];%2006 ; original paper
Pwang2=[2*10^10.48/6.31e11*0.73^2*MvcInMvb,6.31e11/MvcInMvb,1-2.42,1-0.29,1];%vvds
Pwang3=[2*10^10.17/3.21e11*0.73^2*MvcInMvb,3.21e11/MvcInMvb,1-2.42,1-0.29,1];%unified, close to ling
Pwang4=[2*10^10.23/3.43e11*0.73^2*MvcInMvb,3.43e11/MvcInMvb,1-2.56,1-0.34,1];%DR7, 2013
Pguo=[0.129*0.73*MvcInMvb,10^11.4*0.73/MvcInMvb,-0.926,0.261,2.440]; %close to Yang at >M0 ; M200c
Pmoster=[2*0.0282*0.72*MvInMvb,10^11.884*0.72/MvInMvb,-1.057,0.556,1]; %close to ling at low M; tophat Mvir
PmosterScatter=[2*0.02817*0.72*MvInMvb,10^11.899*0.72/MvInMvb,-1.068,0.611,1]; 
Pguofit=[exp(-2.594)*0.73*MvcInMvb,exp(26.32)/MvcInMvb,-1.743,0.592,1]; % refit guo with 4-par model
Pling=[2*10^-1.73,10^11.70,-1.16,0.71,1];%M200b
z=0.;
% PlingZ0=[2*10^-1.73*(1+z)^-0.78,10^(11.70*(1+z)^0.028),-(1.16+0.079*z),0.71*(1+z)^-0.061,1];%M200b
parmosterZ=[0.0351-0.0247*z/(1+z),11.590+1.195*z/(1+z),1.376-0.826*z/(1+z),0.608+0.329*z/(1+z)];
PmosterZ0=[2*parmosterZ(1)*0.72*MvcInMvb,10^parmosterZ(2)*0.72/MvcInMvb,-parmosterZ(3),parmosterZ(4),1];%Moster13, M200c
z=0.2; %redshift evolution makes little difference to the converted Mh(M*) relation.
m=virial_comp(z,0.3);
MvcInMvbZ=m(2)/m(3);
PlingZ=[2*10^-1.73*(1+z)^-0.78,10^(11.70*(1+z)^0.028),-(1.16+0.079*z),0.71*(1+z)^-0.061,1];%M200b
parmosterZ=[0.0351-0.0247*z/(1+z),11.590+1.195*z/(1+z),1.376-0.826*z/(1+z),0.608+0.329*z/(1+z)];
PmosterZ=[2*parmosterZ(1)*0.72*MvcInMvbZ,10^parmosterZ(2)*0.72/MvcInMvbZ,-parmosterZ(3),parmosterZ(4),1];%Moster13, M200c

M=logspace(11,15,100); %for plot
figure;
loglog(M,halo2starP(M,Pling),'r');
hold on;
loglog(M,halo2starP(M,Pwang06),'g');
loglog(M,halo2starP(M,Pwang4),'k')
%%
powerlaw_err=@(x,lnAerr,Berr,corrAB) sqrt(lnAerr^2+(log(x)*Berr).^2+2*log(x)*corrAB*lnAerr*Berr);
powerlawmodel_err=@(m,perr,c) powerlaw_err(m/1e12,perr(1)/log10(exp(1)),perr(2),c); %log(error) for y=10^p(1)*(m/1e12).^p(2)
powerlawmodel=@(m,p) 10^p(1)*(m/1e12).^p(2);
pfiterr=[0.29,0.14];
% pfit=[15.2-0.1,1.64+0.01];
pfit=[15.2,1.64];
cov=0.9;
sm2halo=@(ms) powerlawmodel(ms,pfit); %corrected
%%
par=Pling;
sigma=0.22;
logMmin=7;
logMmax=17;
% sigma0=0.1; %additional scatter due to our SM estimation method
sigma0=0.0;
dndlogMs=@(logms,sigma,par) quadgk(@(logm) halomfunc(logm).*normpdf(logms-log10(halo2starP(10.^logm,par)),0,sqrt(sigma^2+sigma0^2)),logMmin,logMmax);
dndlogMsdlogMh=@(logms,logmh,sigma,par) halomfunc(logmh).*normpdf(logms-log10(halo2starP(10.^logmh,par)),0,sqrt(sigma^2+sigma0^2));
medianMh=@(logMs,sigma,par) quadgk(@(logmh) dndlogMsdlogMh(logMs,logmh,sigma,par)./dndlogMs(logMs,sigma,par).*logmh, logMmin, logMmax);
halo2smfunc=@(mh,par) halomfunc(log10(mh)).*abs(dlnMhdlnMs(mh,par));

ms=logspace(8,12,20);
ywang4=ms;sigma_wang4=0.17; %yang09 and More
ywang=ms;sigma_wang=0.203;
ywang06=ms;sigma_wang06=0.241; %original 2006
ywang2=ms;sigma_wang2=0.203;
ywang3=ms;sigma_wang3=0.240;
yling=ms;sigma_ling=0.22;
ymoster=ms;sigma_moster=0.0; %add stellar dispersion
ymosterZ0=ms;sigma_mosterZ0=0.1;
ymosterScatter=ms;sigma_mosterScatter=0.15;
yguo=ms;sigma_guo=0.0; %add stellar dispersion
ysig=ms;
ynwang=ms;
ynling=ms;
ynmoster=ms;
ynmosterZ0=ms;
figure;
for i=1:numel(ms)
    ywang4(i)=medianMh(log10(ms(i)),sigma_wang4,Pwang4);
    ywang(i)=medianMh(log10(ms(i)),sigma_wang,Pwang);
    ywang06(i)=medianMh(log10(ms(i)),sigma_wang06,Pwang06);
    ywang2(i)=medianMh(log10(ms(i)),sigma_wang2,Pwang2);
    ywang3(i)=medianMh(log10(ms(i)),sigma_wang3,Pwang3);
    yling(i)=medianMh(log10(ms(i)),sigma_ling,Pling);
    ymoster(i)=medianMh(log10(ms(i)),sigma_moster,Pmoster);
    ymosterZ0(i)=medianMh(log10(ms(i)),sigma_mosterZ0,PmosterZ0);
    ymosterScatter(i)=medianMh(log10(ms(i)),sigma_mosterScatter,PmosterScatter);
    yguo(i)=medianMh(log10(ms(i)),sigma_guo,Pguo); 
    ysig(i)=quadgk(@(logmh) dndlogMsdlogMh(log10(ms(i)),logmh,sigma_wang4,Pwang4)./dndlogMs(log10(ms(i)),sigma_wang4,Pwang4).*logmh.^2, logMmin, logMmax);
    ysig(i)=sqrt(ysig(i)-ywang4(i).^2);
    ynwang(i)=dndlogMs(log10(ms(i)),sigma_wang4,Pwang4);
    ynling(i)=dndlogMs(log10(ms(i)),sigma_ling,Pling);
    ynmosterZ0(i)=dndlogMs(log10(ms(i)),sigma_mosterZ0,PmosterZ0);
    ynmoster(i)=dndlogMs(log10(ms(i)),sigma_mosterScatter,PmosterScatter);
    logmh=10:0.1:logMmax;
    plot(logmh,dndlogMsdlogMh(log10(ms(i)),logmh,sigma,par)./dndlogMs(log10(ms(i)),sigma,par));hold on;
end
mh=logspace(12,16);
myfigure;
semilogx(ms,ysig);
xlabel('stellar mass');ylabel('$\sigma_{M_h}$');
myfigure;
loglog(ms,sm2halo(ms),'k','linewidth',5);hold on;
y=powerlawmodel(ms,pfit);
yerrln=powerlawmodel_err(ms,pfiterr,cov);
plot(ms,exp(log(y)+yerrln),'k--');
plot(ms,exp(log(y)-yerrln),'k--');
hold on;
loglog(ms,10.^ywang4,'r');
% loglog(ms,10.^ywang,'r--');
% loglog(ms,10.^ywang06,'r--x');
% loglog(ms,10.^ywang2,'r.-');
% loglog(ms,10.^ywang3,'r:');
loglog(ms,10.^yling,'g');
loglog(ms,10.^ymosterZ0,'b-');
% loglog(ms,10.^ymoster,'b--');
% loglog(ms,10.^ymosterScatter,'b-');
plot(halo2starP(mh,Pguo),mh,'c-');
% loglog(ms,10.^yguo,'c-');
xlabel('stellar mass');ylabel('Median(Mh)');
h=legend('Lensing','err','err','WangLan13','Lingyu13','Moster13','Guo');
set(h,'location','northwest');
% legend('Lensing','err','err','WangLan13','WangLan06','Lingyu13','Moster13','Moster10','Moster10Scatter','Guo');
ylim([1e12,1e16]);
% plot(ms,powerlawmodel(ms,[15.41,1.57]),'y-'); %N>5
% plot(ms,powerlawmodel(ms,[14.74,1.495]),'m-'); %N>1
% plot(ms,powerlawmodel(ms,[15.34,1.47]),'m-');
plot(ms,10.^14.14*(ms/1e12).^(-1.5-1.77*(log10(ms/1e12))),'g-');
% plot(ms,10.^15*(ms/1e12).^(1.0-0.77*(log10(ms/1e12))),'g-');
% print('-depsc','HaloMass-StellarMass.eps')
figure;
% loglog(ms,smfunc(log10(ms)),'kx');hold on;
h1=ploterr(10.^newMF(:,3),newMF(:,4),[],newMF(:,5),'ko','logxy');
hold on;
h2=loglog(ms,ynwang,'r');
h3=loglog(ms,ynling,'g');
h4=loglog(ms,ynmosterZ0,'b');
mh=logspace(7,17);
h5=loglog(halo2starP(mh,Pguo),halo2smfunc(mh,Pguo),'c');
xlabel('stellar mass');ylabel('dN/dlogM');
h=legend([h1(1),h2,h3,h4,h5],'Guo10data','WangLan13','Lingyu13','Moster13','Guo10');
set(h,'location','southwest');
ylim([1e-8,0.1]);
set(gca,'xscale','log');
set(gca,'yscale','log');
% print('-depsc','StellarMassFunc.eps');
%%
% sm2halo=@(ms) 10^15.2*(ms/1e12).^1.64; %raw
sm2halo=@(ms) 10^15.1*(ms/1e12).^1.65; %corrected

dndlogM=@(logmh,sigma) quadgk(@(logm) smfunc(logm).*normpdf(logmh-log10(sm2halo(10.^logm)),0,sigma),0,17);
%%
figure;
mh=logspace(13,15,20);
y00=mh;
y05=mh;
y10=mh;
for i=1:numel(mh)
    y00(i)=dndlogM(log10(mh(i)),0.4);
    y05(i)=dndlogM(log10(mh(i)),0.45);
    y10(i)=dndlogM(log10(mh(i)),0.5);
end
loglog(mh,y00,'r');
hold on;
loglog(mh,y05,'g');
loglog(mh,y10,'b');
hold on;
plot(dndM(:,1),dndM(:,2).*dndM(:,1).*log(10),'k');
%%
figure;
mh=logspace(13,15,20);
y00=mh;
for i=1:numel(mh)
    mmean
end
%%
scale=1;
logh2=0;
% scale=log(10); %convert to dndlnM
% logh2=log10(0.73^2); %convert to Msun
massfun=load('/work/Projects/Lensing/data/StellarMassFunc/massfun.dr72bbright0');
massfun2=load('/work/Projects/Lensing/data/StellarMassFunc/massfun.dr72bbright0.modelmass');
myfigure;
errorbar(massfun(:,3)-logh2,massfun(:,4)/scale,massfun(:,5)/scale,'bo');
hold on;
errorbar(massfun2(:,3)-logh2,massfun2(:,4)/scale,massfun2(:,5)/scale,'rs--');
plot(log10(ms)-logh2,smfunc(log10(ms))/scale,'b-');
plot(log10(ms)-logh2,smfunc(log10(ms),1)/scale,'r-');
set(gca,'yscale','log')
xlim([10,12.3]);
ylim([1e-8,1e-2]);
legend('LW09','Guo10','LW09fit','Guo10fit');
xlabel('log10(Mstar/Msun)');
ylabel('log10(dn/dlnM/[(Mpc/h)$^3$])');
% print('-depsc','stellarmassfun.eps');