cd /work/Projects/Lensing/data/201308
gal=fits_load_bintable('TilingCatv40_kcorr_z00v03.fits');
%%
smdata=fits_load_bintable('StellarMassesv08.fits');
sm_sps=smdata.logmstar+log10(smdata.fluxscale)-2*log10(1.0/0.7);
sm_spsraw=smdata.logmstar-2*log10(1.0/0.7);
i_sps=smdata.absmag_i;
sml_sps=smdata.logmoverl_i;
flag_sps=smdata.fluxscale~=9999&smdata.fluxscale~=-9999&smdata.logmstar~=-99&smdata.logmoverl_i~=-99&smdata.absmag_i~=-99&smdata.nQ>2;
%%
figure;
hist(gal.Z-gal.Z_TONRY)
%%
h=0.7;
dm=@(z) 5*log10(lum_dist(0.3,z)*1e6/h/10);
kcorr=@(p1,p2,p3,p4,p5,z) p1.*z.^4+p2.*z.^3+p3.*z.^2+p4.*z+p5;
zgal=gal.Z_TONRY;
%%
DM=dm(zgal);DM(zgal<=0)=0;
g=gal.G_MODEL-DM-kcorr(gal.PCOEFF_G_1,gal.PCOEFF_G_2,gal.PCOEFF_G_3,gal.PCOEFF_G_4,gal.PCOEFF_G_5,zgal);
i=gal.I_MODEL-DM-kcorr(gal.PCOEFF_I_1,gal.PCOEFF_I_2,gal.PCOEFF_I_3,gal.PCOEFF_I_4,gal.PCOEFF_I_5,zgal);
gevo=g+0.2*zgal;%evolution correction
ievo=i+0.2*zgal;
gal.SM=10.^(1.15+0.70*(g-i)-0.4*i)*0.7^2;%Msun/h^2, what Taylor estimates
gal.SMevo=10.^(1.15+0.70*(gevo-ievo)-0.4*ievo)*0.7^2;%evolution corrected, but is not what Taylor did.
gal.logSML=-0.68+0.7*(g-i)+0.4*4.58;
gal.abs_i=i;
gal.color_gi=g-i;
% clear g i
%%
figure;
plot(sm_sps,sm_spsraw,'.')
%%
gal.logSM_sps=zeros(size(gal.SM));
gal.sml_sps=zeros(size(gal.SM));
gal.i_sps=zeros(size(gal.SM));
gal.flag_sps=zeros(size(gal.SM));
gal.Znew=gal.i_sps;
gal.logSM_spsraw=gal.i_sps;
for i=1:numel(gal.CATAID)
    j=find(smdata.CATAID==gal.CATAID(i));
    if ~isempty(j)
        gal.logSM_sps(i)=sm_sps(j);
        gal.logSM_spsraw(i)=sm_spsraw(j);
        gal.sml_sps(i)=sml_sps(j);
        gal.i_sps(i)=i_sps(j);
        gal.flag_sps(i)=flag_sps(j);
        gal.Znew(i)=smdata.Z_TONRY(j);
    end
end
%%
nbin=20;
myfigure;
h1=plot(smdata.gminusi(flag_sps),sml_sps(flag_sps),'y.','markersize',1);hold on;
f=gal.flag_sps>0&gal.Z_TONRY>0.003;
f=f&gal.NQ>2;
[xx,yy,n,s]=densitygrid(gal.color_gi(f),gal.sml_sps(f),[nbin,nbin],[-0.3,2],[0.8,2.5]);
colormap('cool');
[~,h2]=contour(xx,yy,log10((n+1)),1.2:0.8:4.5);hold on;
% plot(gal.color_gi(f),gal.sml_sps(f),'y.','markersize',1);
hold on;
[xm,ym,ybrd,xmean,ymean,ysig]=skeleton(gal.color_gi(f),gal.sml_sps(f),linspace(-0.5,3,20),0.683);
h3=errorbar(xm,ymean,ysig,'o');
h4=plot(xm,ym,'r-');
plot(xm,ybrd(:,1),'r-');
plot(xm,ybrd(:,2),'r-');
h5=plot(xm,(-0.68+0.7*xm+0.4*4.58),'k-','linewidth',4);
xlim([-0.15,1.7])
xlabel('g-i');ylabel('log(M/L)-log(M$_\odot/L_{AB}$)');
h=legend([h1,h2,h3,h4,h5],'SPS output','SDSS color','SDSS mean','SDSS percentiles','Fit Relation');
set(h,'location','southeast');
print('-depsc','/work/Projects/Lensing/outputv4/SM-color.eps');
%%
f=gal.flag_sps>0&gal.Z_TONRY>0.003;
f=f&gal.NQ>2;
figure;
plot(gal.logSML(f),gal.sml_sps(f),'y.','markersize',1);
hold on;
hold on;
[xm,ym,ybrd,xmean,ymean,ysig]=skeleton(gal.logSML(f),gal.sml_sps(f),linspace(0.5,2.5,20),0.683);
h3=errorbar(xm,ymean,ysig,'o');
h4=plot(xm,ym,'r-');
plot(xm,ybrd(:,1),'r-');
plot(xm,ybrd(:,2),'r-');
xlim([0.5,2.5]);ylim([0.5,2.5]);
plot([0.5,2.5],[0.5,2.5],'k-','linewidth',3);
%%
myfigure;
% logSM=log10(gal.SM);
zbins=[0.01,0.14;0.14,0.25;0.25,0.35;0.35,0.55];
% cbins=[0.,0.4;0.4,0.8;0.8,1.2;1.2,1.6];
for i=1:4
subplot(2,2,i);
f0=gal.Z_TONRY>zbins(i,1)&gal.Z_TONRY<zbins(i,2);
% f0=gal.color_gi>cbins(i,1)&gal.Z_TONRY<cbins(i,2);
x=log10(gal.SM(f&f0));y=gal.logSM_sps(f&f0)-x;
% x=gal.logSM_sps(f&f0);y=gal.logSM_sps(f&f0)-log10(gal.SM(f&f0));
plot(x,y,'c.','markersize',1);hold on;
[xm,ym,ybrd,xmean,ymean,ysig]=skeleton(x,y,linspace(8.5,11.5,10),0.683);
hold on;
h3=errorbar(xm,ymean,ysig,'o');
h4=plot(xm,ym,'r-');
plot(xm,ybrd(:,1),'r-');
plot(xm,ybrd(:,2),'r-');
% plot([6,14],[6,14],'k-');
plot([6,14],[0,0],'k--');
plot([6,14],[0.1,0.1],'k--');
xlim([8.5,11.5]);
ylim([-0.5,0.5]);
title([num2str(zbins(i,1)),'$<z<$',num2str(zbins(i,2))]);
xlabel('log(Mcolor)');ylabel('log(Msps/Mcolor)');
end
% print('-depsc','/work/Projects/Lensing/outputv4/SMRat-SMcolor-z.eps');
%%
myfigure;
logSM=log10(gal.SM);
zbins=[0.01,0.14;0.14,0.25;0.25,0.35;0.35,0.55];
cbins=[0.,0.4;0.4,0.8;0.8,1.2;1.2,1.6];
for i=1:4
subplot(2,2,i);
f0=gal.Z_TONRY>zbins(i,1)&gal.Z_TONRY<zbins(i,2);
% f0=gal.color_gi>cbins(i,1)&gal.Z_TONRY<cbins(i,2);
plot(logSM(f&f0),gal.logSM_sps(f&f0),'.','markersize',1);hold on;
[xm,ym]=skeleton(logSM(f&f0),gal.logSM_sps(f&f0),linspace(8,12,10),0.683);
% plot(mean(logSM(f&f0)),mean(gal.logSM_sps(f&f0)),'rx')
10.^(mean(logSM(f&f0))-mean(gal.logSM_sps(f&f0)))
hold on;
plot(xm,ym,'r-');
plot([6,14],[6,14],'k-');
xlim([7,13]);
ylim([7,13]);
title([num2str(zbins(i,1)),'$<z<$',num2str(zbins(i,2))]);
xlabel('log(Mcolor)');ylabel('log(Msps)');
end
print('-depsc','/work/Projects/Lensing/outputv4/SMcolor-SMsed.eps');
myfigure;
logSM=log10(gal.SMevo/0.7^2);
zbins=[0.01,0.14;0.14,0.25;0.25,0.35;0.35,0.55];
for i=1:4
subplot(2,2,i);
f0=gal.Z_TONRY>zbins(i,1)&gal.Z_TONRY<zbins(i,2);
plot(logSM(f&f0),gal.logSM_sps(f&f0),'.','markersize',1);
hold on;
[xm,ym]=skeleton(logSM(f&f0),gal.logSM_sps(f&f0),linspace(8,12,20),0.683);
plot(xm,ym,'r-');
plot([6,14],[6,14],'k-');
xlim([7,13]);
ylim([7,13]);
title([num2str(zbins(i,1)),'$<z<$',num2str(zbins(i,2))]);
xlabel('log(McolorEvo)');ylabel('log(Msps)');
end
print('-depsc','/work/Projects/Lensing/outputv4/SMcolorEvo-SMsed.eps');
%%
f=gal.flag_sps>0&gal.Z_TONRY>0.003;
f=f&gal.NQ>2;
myfigure; 
nbin=50;
% x=log10(gal.SM(f));y=gal.logSM_sps(f)-x;
x=gal.logSM_sps(f);y=x-log10(gal.SM(f));
plot(x,y,'c.','markersize',1);hold on;
% [xx,yy,n,s]=densitygrid(x,y,[nbin,nbin],[7,13],[0.5,1.5]);
% colormap('cool');
% [~,h2]=contour(xx,yy,log10((n+1)),0.6:0.5:5);hold on;
[xm,ym,ybrd,xmean,ymean,ysig]=skeleton(x,y,linspace(7,12,15),0.683);
hold on;
h3=errorbar(xm,ymean,ysig,'o');
h4=plot(xm,ym,'r-');
plot(xm,ybrd(:,1),'r-');
plot(xm,ybrd(:,2),'r-');
% ploterr(xm,ym,[],{ybrd(:,1),ybrd(:,2)},'ro');
% plot([6,14],[6,14],'k-');
plot([6,14],[0,0],'k--');
plot([6,14],[0.1,0.1],'k--');
xlim([7,13]);
ylim([-0.5,0.5]);
xlabel('log(Msps)');
ylabel('log(Msps)-log(Mcolor)');
set(gca,'ytick',-0.5:0.1:0.5);
print('-depsc','/work/Projects/Lensing/outputv4/SMRat-SMsps.eps');
%%
figure;
f=gal.flag_sps>0&gal.Z_TONRY>0.003;
f=f&gal.NQ>2;
f=f&logSM>9&logSM<10;
linhist(gal.logSM_sps(f)-logSM(f),-2:0.05:2,'stairsnorm');
hold on;
[m,s]=normfit(gal.logSM_sps(f)-logSM(f));
plot(-2:0.05:2,normpdf(-2:0.05:2,m,s)*0.05,'r');
%%
figure;
% f0=abs(gal.Znew-gal.Z_TONRY)<0.001;
f0=gal.NQ>0;
plot(gal.i_sps(f&f0),gal.abs_i(f&f0),'.','markersize',1);hold on;
% xlim([-30,-5]);
% axis equal
plot([-30,-5],[-30,-5],'k');
%%
[xx,yy,n,s]=densitygrid(gal.i_sps(f&f0),gal.abs_i(f&f0),[nbin,nbin]);
colormap('cool');
[~,h2]=contour(xx,yy,log10((n+1)));hold on;
%%
figure;
plot(gal.Znew(f),gal.Z_TONRY(f),'.');
%%
figure;
loghist(gal.Znew(f)./gal.Z_TONRY(f),logspace(-3,3,10),'norm')
%%
MvInMvb=0.893;
MvcInMvb=0.733;
halo2star=@(M,A,M0,alpha,beta,gamma) A./((M/M0).^alpha+(M/M0).^beta).^gamma.*M;
halo2starP=@(M,P) halo2star(M,P(1),P(2),P(3),P(4),P(5)); %Msun/h^2 for Mstar, Msun/h for Mhalo
Pwang=[2*10^10.35/4e11*0.73^2*MvcInMvb,4e11/MvcInMvb,1-2.42,1-0.29,1];%2006 ; probably M200c, convert to M200b
Pwang2=[2*10^10.48/6.31e11*0.73^2*MvcInMvb,6.31e11/MvcInMvb,1-2.42,1-0.29,1];%vvds
Pwang3=[2*10^10.17/3.21e11*0.73^2*MvcInMvb,3.21e11/MvcInMvb,1-2.42,1-0.29,1];%unified, close to ling
Pwang4=[2*10^10.23/3.43e11*0.73^2*MvcInMvb,3.43e11/MvcInMvb,1-2.56,1-0.34,1];%DR7, 2013
Pguo=[0.129*0.73*MvcInMvb,10^11.4*0.73/MvcInMvb,-0.926,0.261,2.440]; %close to Yang at >M0 ; M200c
Pmoster=[2*0.0282*0.72*MvInMvb,10^11.884*0.72/MvInMvb,-1.057,0.556,1]; %close to ling at low M; tophat Mvir
Pguofit=[exp(-2.594)*0.73*MvcInMvb,exp(26.32)/MvcInMvb,-1.743,0.592,1]; % refit guo with 4-par model
Pling=[2*10^-1.73,10^11.70,-1.16,0.71,1];%M200b
Pfit=[10^8.4 10^9.6 -3.2 -0.6];
Pfit=[Pfit(1)*2/Pfit(2),Pfit(2),1+Pfit(3),1+Pfit(4),1];%M_*1=A*M_1/2, alpha-1, beta-1
% Pyang=[2*0.0384*0.72,10^12.067*0.72,-0.71,0.698,1]; %fitted by Moster;
% M180b. deviate from the original. don't trust it.
z=0.;
% PlingZ0=[2*10^-1.73*(1+z)^-0.78,10^(11.70*(1+z)^0.028),-(1.16+0.079*z),0.71*(1+z)^-0.061,1];%M200b
parmosterZ=[0.0351-0.0247*z/(1+z),11.590+1.195*z/(1+z),1.376-0.826*z/(1+z),0.608+0.329*z/(1+z)];
PmosterZ0=[2*parmosterZ(1)*0.72*MvcInMvb,10^parmosterZ(2)*0.72/MvcInMvb,-parmosterZ(3),parmosterZ(4),1];%Moster13, M200c

z=0.2;
m=virial_comp(z,0.3);
MvcInMvbZ=m(2)/m(3);
PlingZ=[2*10^-1.73*(1+z)^-0.78,10^(11.70*(1+z)^0.028),-(1.16+0.079*z),0.71*(1+z)^-0.061,1];%M200b
parmosterZ=[0.0351-0.0247*z/(1+z),11.590+1.195*z/(1+z),1.376-0.826*z/(1+z),0.608+0.329*z/(1+z)];
PmosterZ=[2*parmosterZ(1)*0.72*MvcInMvbZ,10^parmosterZ(2)*0.72/MvcInMvbZ,-parmosterZ(3),parmosterZ(4),1];%Moster13, M200c

M=logspace(11,15,100); %for plot
Mswang=halo2starP(M,Pwang);
Msguo=halo2starP(M,Pguo);
Msmoster=halo2starP(M,Pmoster);
Msling=halo2starP(M,Pling);
Msyang=10^10.86*(M/10^12.08).^(0.22+1.61)./(1+M/10^12.08).^1.61; %M180b

myfigure;
loglog(M,Msguo,'g',M,Msmoster,'b',M,Msling,'k');
hold on;
% loglog(M,Mswang,'r:');
% loglog(M,halo2starP(M,Pwang2),'r--');
% loglog(M,halo2starP(M,Pwang3),'r-');
loglog(M,halo2starP(M,Pwang4),'r-');
% loglog(M,halo2starP(M,Pguofit),'g--');
% loglog(M,halo2starP(M,Pyang),'c');
loglog(M,Msyang,'c--')
loglog(M,halo2starP(M,PlingZ),'k--');
loglog(M,halo2starP(M,PmosterZ),'b--');
loglog(M,halo2starP(M,PmosterZ0),'b-');
% loglog(M,halo2starP(M,Pfit),'ko');
ms=logspace(10,11.5);
% p5=[15.2-0.1,1.64+0.01];
plot(10^15.1*(ms/1e12).^1.65,ms,'rd');
% f=grp.Mult>2;
% plot(starM(f),grp.IterCenSM(f),'r.')
l=legend('Guo10','Moster10','Wang13','WangLan13','Yang08');
set(l,'location','southeast');
xlabel('$M_{200b} [M_\odot/h]$');
ylabel('$M_\star [M_\odot/h^2]$');
% print('-depsc','/work/Projects/Lensing/outputv4/StellarMass-HaloMass.eps');


figure;
loglog(M,Msguo./M,'g',M,Msmoster./M,'b',M,Msling./M,'k');
hold on;
loglog(M,halo2starP(M,Pwang4)./M,'r');
loglog(M,Msyang./M,'c')
% figure;
% semilogx(M,halo2starP(M,Pguofit)./Msguo-1)
%% try to invert: function form ok, but pars differ
newp=@(p) [p(1)*p(2),p(2),1-p(3),1-p(4),1];

invp=@(p) [p(2)/2,p(1)/2,1/p(3),1/p(4)];

star2haloP=@(M,P) P(1)*((M/P(2)).^P(3)+(M/P(2)).^P(4));
minv=halo2starP(M,Pling);
loglog(M,minv,'r');
hold on;
plot(star2haloP(minv,invp(newp(Pling))),minv,'g')
%%
grp=fits_load_bintable('G3Cv6/G3CFoFGroupv06.fits');
grp.IterCenSM=zeros(size(grp.GroupID));
grp.BCGSM=grp.IterCenSM;
for i=1:numel(grp.GroupID)
    grp.IterCenSM(i)=gal.SM(gal.CATAID==grp.IterCenCATAID(i));
    grp.BCGSM(i)=gal.SM(gal.CATAID==grp.BCGCATAID(i));
end
% SMall=zeros(max(gal.CATAID));
% SMall(gal.CATAID)=gal.SM;
% grp.IterCenSM=SMall(grp.IterCenCATAID);
% grp.BCGSM=SMall(grp.BCGCATAID);
%%
M=logspace(9,17,100); %for interpolation
Msyang=10^10.86*(M/10^12.08).^(0.22+1.61)./(1+M/10^12.08).^1.61; %M180b
grp.HODMassIterMoster=spline(halo2starP(M,Pmoster),M,grp.IterCenSM);
grp.HODMassBCG=spline(halo2starP(M,Pmoster),M,grp.BCGSM);
grp.HODMassIterGuo=spline(halo2starP(M,Pguo),M,grp.IterCenSM);
grp.HODMassIterWang=spline(halo2starP(M,Pwang),M,grp.IterCenSM);
grp.HODMassIterLing=spline(halo2starP(M,Pling),M,grp.IterCenSM);
grp.HODMassIterYang=spline(Msyang,M,grp.IterCenSM);
%compare M_HOD with M_dyn/lum
%%
bins=logspace(9,17,10);
figure;
loglog(grp.MassAfunc,grp.HODMassIter,'y.')
hold on;
[xm,ym,ylim]=skeleton(grp.MassAfunc,grp.HODMassIter,bins,0.683);
plot(xm,ym,'r');
plot(xm,ylim(:,1),'r--');
plot(xm,ylim(:,2),'r--');
% loglog(grp.MassAfunc,grp.HODMassBCG,'g.')
% [xm,ym,ylim]=skeleton(grp.MassAfunc,grp.HODMassBCG,bins,0.683);
% plot(xm,ym,'g');
% plot(xm,ylim(:,1),'g--');
% plot(xm,ylim(:,2),'g--');
% loglog(grp.MassAfunc,grp.HODMassIterGuo,'b.')
[xm,ym,ylim]=skeleton(grp.MassAfunc,grp.HODMassIterGuo,bins,0.683);
plot(xm,ym,'b');
plot(xm,ylim(:,1),'b--');
plot(xm,ylim(:,2),'b--');
% loglog(grp.MassAfunc,grp.HODMassIterWang,'c.')
[xm,ym,ylim]=skeleton(grp.MassAfunc,grp.HODMassIterWang,bins,0.683);
plot(xm,ym,'c');
plot(xm,ylim(:,1),'c--');
plot(xm,ylim(:,2),'c--');
% loglog(grp.MassAfunc,grp.MassA,'c.')
plot([1e9,1e15],[1e9,1e15],'k')
%%
bins=logspace(10,16,15);
f=grp.Mult>2;
myfigure;
% loghist(grp.MassAfunc(f),bins,'plot','r')
loghist(grp.LumMass(f),bins,'plot','r-o')
hold on;
loghist(grp.MassProxy(f),bins,'plot','r-s')
loghist(grp.HODMassIterMoster(f),bins,'plot','g')
loghist(grp.HODMassIterGuo(f),bins,'plot','b')
loghist(grp.HODMassIterWang(f),bins,'plot','c')
loghist(grp.HODMassIterLing(f),bins,'plot','k')
% loghist(grp.HODMassBCG(f),bins,'plot','g--')
loghist(grp.HODMassIterYang(f),bins,'plot','m')
% loghist(spline(halo2starP(M,Pguofit),M,grp.IterCenSM(f)),bins,'plot','b--')
legend('lum','dyn','moster','guo','wang10u','wang13','yang')
xlabel('$M_{200b} [M_\odot/h]$');
ylabel('counts');
print('-depsc','/work/Projects/Lensing/outputv4/HODMass.eps');