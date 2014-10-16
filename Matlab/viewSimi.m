%% to compare similar quantities
%%
load_G3Cv4;
%%
nbin=20;
%% Mass vs. LumMass
nbin=20;
f=grpmock.Mult>1;
myfigure;
[x,y,dyx]=loghist(grpmock.Mass(f)./grpmock.LumMass(f),nbin,'stairs','--');
% loglog(x,dyx,'--');
hold on;

f=grp.Mult>1;
[x,y,dyx]=loghist(grp.Mass(f)./grp.LumMass(f),nbin,'stairs','-');
% loglog(x,dyx,'-');

l=legend('Mock','GAMA');
set(l,'interpreter','latex','location','northwest','fontsize',15);
% set(gca,'xlim',[0,1]);
xlabel('$M_{dyn}/M_{lum}$');
ylabel('$dN$');
print('-depsc','GroupMasses.eps');

myfigure;
f=grpmock.Mult>1;
loglog(grpmock.Mass(f),grpmock.LumMass(f),'o','markersize',2);
hold on;
f=grp.Mult>1;
loglog(grp.Mass(f),grp.LumMass(f),'r.','markersize',2);
plot([1e10,1e16],[1e10,1e16],'r-');
xlabel('$M_{dyn}/(M_{\odot}/h)$');
ylabel('$M_{lum}/(M_{\odot}/h)$');
l=legend('Mock','GAMA','y=x');set(l,'location','northwest','fontsize',15);
set(gca,'xlim',[1e10,1e16],'ylim',[1e10,1e16]);
print('-depsc','GroupMassLumMass.eps');
%% SMInt vs. Mass
myfigure;
loglog(grp.TotSMInt,grp.Mass,'.','markersize',4);hold on;plot([1e8,1e14],[1e8,1e14],'r-');
xlabel('Extrapolated $M_{star}/(M_{\odot}/h)$');ylabel('$M_{dyn}/(M_{\odot}/h)$');
l=legend('GAMA','y=x');set(l,'location','southeast');
set(gca,'xlim',[1e8,1e16],'ylim',[1e8,1e16]);
print('-depsc','GroupStarDynmMasses.eps');
%% SMInt vs. LumMass
myfigure;
loglog(grp.TotSMInt,grp.LumMass,'.','markersize',4);hold on;plot([1e7,1e14],[1e7,1e14],'r-');
xlabel('Extrapolated $M_{star}/(M_{\odot}/h)$');ylabel('$M_{lum}/(M_{\odot}/h)$');
l=legend('GAMA','y=x');set(l,'location','northwest');
set(gca,'xlim',[1e8,1e16],'ylim',[1e8,1e16]);
print('-depsc','GroupStarLumMasses.eps');
%% SM vs. SMInt
myfigure;
loglog(grp.TotSM,grp.TotSMInt,'.','markersize',4);hold on;plot([1e7,1e14],[1e7,1e14],'r-');
xlabel('$M_{star}/(M_{\odot}/h)$');ylabel('Extrapolated $M_{star}/(M_{\odot}/h)$');
l=legend('GAMA','y=x');set(l,'location','northwest');
print('-depsc','GroupStars.eps');
%% Flux vs. FluxInt
myfigure;
plot(grp.TotFlux,grp.TotFluxInt,'.','markersize',4);hold on;plot([-26,-14],[-26,-14],'r-');
xlabel('$M_{r}$');ylabel('Extrapolated $M_{r}$');
l=legend('GAMA','y=x');set(l,'location','northwest');
print('-depsc','GroupMags.eps');
%% Rad50, RadSig, vs. Rad100
nbin=25;
f=grpmock.Mult>5;
myfigure;
linhist(grpmock.Rad50(f)./grpmock.Rad1Sig(f),nbin,'stairs','r--');
hold on;
linhist(grpmock.Rad50(f)./grpmock.Rad100(f),nbin,'stairs','g--');
f=grp.Mult>5;
linhist(grp.Rad50(f)./grp.Rad1Sig(f),nbin,'stairs','r-');
hold on;
linhist(grp.Rad50(f)./grp.Rad100(f),nbin,'stairs','g-');
l=legend('Mock $R_{\sigma}$','Mock $R_{max}$','GAMA $R_{\sigma}$','GAMA $R_{max}$');
set(l,'interpreter','latex','location','northwest','fontsize',15);
set(gca,'xlim',[0,1]);
xlabel('$R_{0.5}/R_{?}$');
ylabel('$dN$');
title('$Multiplicity>5$');

print('-depsc','GroupSizeRat.eps');

myfigure;
f=grp.Mult>5;
plot(grp.Rad1Sig(f),grp.Rad50(f),'.');
hold on;
plot(grp.Rad1Sig(f),grp.Rad100(f),'g.');
plot([0,1.],[0,1.],'r-');
xlabel('$R_{\sigma}$');ylabel('$R_{.50}$ or $R_{max}$');
title('GAMA $Multiplicity>5$');
print('-depsc','GroupSizes.eps');
%% centers
datadir='/home/kam/Projects/Lensing/data/';
% global cat
% grp.HybMass=grp.LumMass;
% [cat{1},cat{2},cat{3}]=split_mockcat(grp);
% % labels=mark_mz_bin('halomass');
% grp=structcat([cat{1},cat{2},cat{3}]);

MultMin=1;
f=grp.Mult>MultMin;
figure;
plot(grp.IterCenRef(f),grp.BCGRef(f),'.');
DCenBCG=acos(ccdist([grp.BCGRA,grp.BCGDEC],[grp.CenRA,grp.CenDEC])).*AD_dist_flat(0.3,0,grp.Zfof)./grp.Rad1Sig;
% DCenIter=acos(ccdist([grp.IterCenRA,grp.IterCenDEC],[grp.CenRA,grp.CenDEC])).*AD_dist_flat(0.3,0,grp.Zfof)./grp.Rad1Sig;
DCenIter=acos(ccdist([grp.IterCenRA,grp.IterCenDEC],[grp.CenRA,grp.CenDEC])).*AD_dist_flat(0.3,0,grp.Zfof)./grp.Rad50;
DBCGIter=acos(ccdist([grp.BCGRA,grp.BCGDEC],[grp.IterCenRA,grp.IterCenDEC])).*AD_dist_flat(0.3,0,grp.Zfof)./grp.Rad50;

myfigure;
x=0:0.02:1;
linhist(DCenBCG(f),x,'stairsnorm','r');hold on;
linhist(DCenIter(f),x,'stairsnorm','g');
[xm,ym]=linhist(DBCGIter(f),x,'stairsnorm','k');
legend('Cen-BCG','Cen-Iter','BCG-Iter');
xlabel('$\Delta R/R_{\sigma}$');ylabel('Probability');title(['$GAMA Mult>',num2str(MultMin),'$']);
% print('-depsc','/home/kam/Projects/Lensing/outputv4/GAMACenters.eps');
%%
DCenIter=acos(ccdist([grp.IterCenRA,grp.IterCenDEC],[grp.CenRA,grp.CenDEC])).*AD_dist_flat(0.3,0,grp.Zfof)./grp.Rad50;%.*(1+grp.Zfof);
x=linspace(0,1,15);
myfigure;
m=[0,5e12;5e12,3.3e13;3.3e13,1e16];
mlt=[0,3;3,5;5,10;10,1000];
z=[0,0.1;0.1,0.2;0.2,0.3;0.3,0.5];
h=[];s=[];
colors='krgbcm';
symb='*os>';
for i=2:4
%     ff=grp.LumMass>m(i,1)&grp.LumMass<m(i,2);
    ff=grp.Mult>=mlt(i,1)&grp.Mult<mlt(i,2)&DCenIter<1;
% ff=grp.Mult>2&grp.Zfof>z(i,1)&grp.Zfof<z(i,2);
%     ff=ff&f;
% ff=f;
    [xm,ym,dyx]=linhist(DCenIter(ff),x);%hold on;
    htmp=plot(xm,dyx./sum(ym),symb(i),'color',colors(i),'markersize',10);hold on;
    h=[h,htmp];
    sig=raylfit(DCenIter(ff),0);
    s=[s,sig];
%     htmp=plot(xm,raylpdf(xm,sig),'color',colors(i));
%     h=[h,htmp];
% loglog(xm,ym,'.-','color',colors(i));hold on;
end
sig=raylfit(DCenIter(grp.Mult>2&DCenIter<1),0);
htmp=plot(x,raylpdf(x,sig),'k');
h=[h,htmp];
xlabel('$D_{CoL-Iter}/R_{50}$');ylabel('Probability Density');
l=legend(h,'$N=3\sim4$','$N=5\sim9$','$N>9$','$N>2$ fit');set(l,'interpreter','latex');
% gtext(num2str(s,'%.2f\n'))
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/DisplacementDistribution.eps');
%%
theta=grp.BCGRA;
fai=theta;
for i=1:numel(grp.BCGRA)
    [theta(i),fai(i)]=sphshift2(grp.BCGRA(i),grp.BCGDEC(i),grp.CenRA(i),grp.CenDEC(i));
end
Y=[theta.*AD_dist_flat(0.3,0,grp.Zfof)./grp.Rad1Sig,fai.*AD_dist_flat(0.3,0,grp.Zfof)./grp.Rad1Sig]*pi/180;
% fai=atan2(grp.BCGDEC-grp.CenDEC,cosd(grp.CenDEC).*(grp.BCGRA-grp.CenRA));
% X=[DCenBCG.*cos(fai),DCenBCG.*sin(fai)];
% obj=gmdistribution.fit(X,2);
[mu,sig]=normfit(Y);
%%
figure;
histfit(Y(:,1),50)
%%
myfigure;
% x=0:0.1:2.5;
x=logspace(-2,0,13);
colors='rgbk';
for i=2:4
ff=grp.Mult==i;
ff=ff;
[xm,ym,dyx]=loghist(DCenBCG(ff),x);%hold on;
semilogx(xm,dyx./sum(ym),'.-','color',colors(i));hold on;
% loglog(xm,ym,'.-','color',colors(i));hold on;
end
xlabel('comoving BCG-CoL Displacement (Mpc/h)');ylabel('Fraction');
legend('1','2','3','4');
% print('-depsc','/home/kam/Projects/Lensing/outputv4/GAMACenters_bin.eps');

%%
f=grpmock.Mult>MultMin;
figure;
plot(grpmock.IterCenRef(f),grpmock.BCGRef(f),'.');
DCenBCG=acos(ccdist([grpmock.BCGRA(f),grpmock.BCGDEC(f)],[grpmock.CenRA(f),grpmock.CenDEC(f)])).*AD_dist_flat(0.3,0,grpmock.Zfof(f))./grpmock.Rad1Sig(f);
DCenIter=acos(ccdist([grpmock.IterCenRA(f),grpmock.IterCenDEC(f)],[grpmock.CenRA(f),grpmock.CenDEC(f)])).*AD_dist_flat(0.3,0,grpmock.Zfof(f))./grpmock.Rad1Sig(f);
DBCGIter=acos(ccdist([grpmock.BCGRA(f),grpmock.BCGDEC(f)],[grpmock.IterCenRA(f),grpmock.IterCenDEC(f)])).*AD_dist_flat(0.3,0,grpmock.Zfof(f))./grpmock.Rad1Sig(f);
myfigure;
x=0:0.1:2.5;
linhist(DCenBCG,x,'stairsnorm','r');hold on;
linhist(DCenIter,x,'stairsnorm','g');
linhist(DBCGIter,x,'stairsnorm','k');
legend('Cen-BCG','Cen-Iter','BCG-Iter');
xlabel('$\Delta R/R_{\sigma}$');ylabel('Probability');title(['$Mock Mult>',num2str(MultMin),'$']);
    