Continuum SN cut is not the same as EMI_SN cut
HalphaFlux cut is more close to EMI_SN cut
Or perhaps combine EMI_SN with Cont_SN
%%
xbin=logspace(-2,5);
loghist(gal.EmLineSN(gal.FlagSFRFluxCut),xbin,'-','ro')
hold on;
loghist(gal.EmLineSN(gal.FlagSFRContCut),xbin,'-','gd')
loghist(gal.EmLineSN(gal.FlagSFRCombCut),xbin,'-','bs')
loghist(gal.EmLineSN(gal.IsAGN~=1&gal.Zspec<0.31&gal.Zspec>-0.01&(gal.FlagWarn==4|gal.FlagWarn==-1)),xbin,'-','k')
legend('FluxCut','ContCut','CombCut','Basic')
xlabel('EMI_SN')
%%
xbin=logspace(-2,5);
loghist(gal.ContSN(gal.FlagSFRFluxCut),xbin,'-','ro')
hold on;
loghist(gal.ContSN(gal.FlagSFRContCut),xbin,'-','gd')
loghist(gal.ContSN(gal.FlagSFRCombCut),xbin,'-','bs')
loghist(gal.ContSN(gal.IsAGN~=1&gal.Zspec<0.31&gal.Zspec>-0.01&(gal.FlagWarn==4|gal.FlagWarn==-1)),xbin,'-','k')
legend('FluxCut','ContCut','CombCut','Basic')
xlabel('CONT_SN')
%%
xbin=logspace(-2,5);
loghist(gal.FluxHalpha(gal.FlagSFRFluxCut),xbin,'-','ro')
hold on;
loghist(gal.FluxHalpha(gal.FlagSFRContCut),xbin,'-','gd')
loghist(gal.FluxHalpha(gal.FlagSFRCombCut),xbin,'-','bs')
loghist(gal.FluxHalpha(gal.IsAGN~=1&gal.Zspec<0.31&gal.Zspec>-0.01&(gal.FlagWarn==4|gal.FlagWarn==-1)),xbin,'-','k')
legend('FluxCut','ContCut','CombCut','Basic')
xlabel('FluxHalpha')
%%
smth=@(n) conv2(n,[0.02,0.05,0.02;0.05,0.72,0.05;0.02,0.05,0.02;],'same');
nbin=50;
percents=[0.3,0.6,0.9];
figure;
f=gal.IsAGN~=1&gal.Zspec<0.31&gal.Zspec>-0.01&gal.EmLineSN>0&gal.ContSN>0;%&(gal.FlagWarn==4|gal.FlagWarn==-1);
x=log10(gal.EmLineSN(f));
y=log10(gal.ContSN(f));
plot(x,y,'.','markersize',1); hold on;
[xx,yy,n,s]=densitygrid(x,y,[nbin,nbin],[-2,5],[-2,2]);
[~,h2]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'k');hold on;
%%
smth=@(n) conv2(n,[0.02,0.05,0.02;0.05,0.72,0.05;0.02,0.05,0.02;],'same');
nbin=50;
percents=[0.3,0.6,0.9];
figure;
f=gal.IsAGN~=1&gal.Zspec<0.31&gal.Zspec>-0.01&gal.EmLineSN>0;%&(gal.FlagWarn==4|gal.FlagWarn==-1);
x=log10(gal.EmLineSN(f));
y=log10(gal.FluxHalpha(f));
plot(x,y,'.','markersize',1); hold on;
[xx,yy,n,s]=densitygrid(x,y,[nbin,nbin],[-2,5],[0,5]);
[~,h2]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'k');hold on;
%%
smth=@(n) conv2(n,[0.02,0.05,0.02;0.05,0.72,0.05;0.02,0.05,0.02;],'same');
nbin=50;
percents=[0.3,0.6,0.9];
figure;
f=gal.IsAGN~=1&gal.Zspec<0.31&gal.Zspec>-0.01&gal.EmLineSN>0;%&(gal.FlagWarn==4|gal.FlagWarn==-1);
x=log10(gal.EmLineSN(f));
y=log10(gal.LumHalpha(f));
plot(x,y,'.','markersize',1); hold on;
[xx,yy,n,s]=densitygrid(x,y,[nbin,nbin],[-2,5],[30,40]);
[~,h2]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'k');hold on;
%%
% smth=@(n) conv2(n,[0.05,0.1,0.05;0.1,0.4,0.1;0.05,0.1,0.05;],'same');
smth=@(n) conv2(n,[0.02,0.05,0.02;0.05,0.72,0.05;0.02,0.05,0.02;],'same');
% smth=@(n) n;
nbin=50;
percents=[0.3,0.6,0.9];
figure;
f=gal.IsAGN~=1&gal.Zspec<0.31&gal.Zspec>-0.01&gal.EmLineSN>0;%&(gal.FlagWarn==4|gal.FlagWarn==-1);
x=log10(gal.EmLineSN(f));
y=log10(gal.SSFR(f));
plot(x,y,'.','markersize',1); hold on;
plot([-5,5],[-1.5,-1.5],'r-')
[xx,yy,n,s]=densitygrid(x,y,[nbin,nbin],[-2,5],[-3,1]);
[~,h2]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'k');hold on;
xlim([-2,4]);ylim([-3,1]);
xlabel('log(LineSN)');ylabel('log(SSFR)');
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/LineSN-SSFR.eps')
%%
% smth=@(n) conv2(n,[0.05,0.1,0.05;0.1,0.4,0.1;0.05,0.1,0.05;],'same');
smth=@(n) conv2(n,[0.02,0.05,0.02;0.05,0.72,0.05;0.02,0.05,0.02;],'same');
% smth=@(n) n;
nbin=50;
percents=[0.3,0.6,0.9];
figure;
f=gal.IsAGN~=1&gal.Zspec<0.31&gal.Zspec>-0.01&gal.EmLineSN>0;%&(gal.FlagWarn==4|gal.FlagWarn==-1);
x=log10(gal.FluxHalpha(f));
y=log10(gal.SSFR(f));
plot(x,y,'.','markersize',1); hold on;
plot([-1,3],[-1.5,-1.5],'r-')
[xx,yy,n,s]=densitygrid(x,y,[nbin,nbin],[-1,5],[-5,5]);
[~,h2]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'k');hold on;
xlim([0.5,4]);ylim([-4,3]);
xlabel('FluxHalpha');ylabel('SSFR');
%%
% smth=@(n) conv2(n,[0.05,0.1,0.05;0.1,0.4,0.1;0.05,0.1,0.05;],'same');
smth=@(n) conv2(n,[0.02,0.05,0.02;0.05,0.72,0.05;0.02,0.05,0.02;],'same');
% smth=@(n) n;
nbin=50;
percents=[0.3,0.6,0.9];
figure;
f=gal.IsAGN~=1&gal.Zspec<0.31&gal.Zspec>-0.01&(gal.FlagWarn==4|gal.FlagWarn==-1);
x=log10(gal.LumHalpha(f));
y=log10(gal.SSFR(f));
plot(x,y,'.','markersize',1); hold on;
plot([30,40],[-1.5,-1.5],'r-')
[xx,yy,n,s]=densitygrid(x,y,[nbin,nbin],[30,40],[-5,5]);
[~,h2]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'k');hold on;
xlim([31,36]);ylim([-3,1]);
xlabel('LumHalpha');ylabel('SSFR');
%%
smth=@(n) conv2(n,[0.02,0.05,0.02;0.05,0.72,0.05;0.02,0.05,0.02;],'same');
nbin=30;
percents=[0.3,0.6,0.9];
figure;
f=gal.FlagSFR>=2&gal.CentralSampleIter>0;
x=log10(gal.SMsps(f));
y=(gal.Zspec(f));
% plot(x,y,'.','markersize',1); hold on;
[xx,yy,n,s]=densitygrid(x,y,[nbin,nbin],[8,13],[0,0.4]);
[~,h2]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'k');hold on;
f=gal.FlagSFR>=2&gal.CentralSampleIter>0;
f=gal.FlagSFRBaseCut>0&~f;
x=log10(gal.SMsps(f));
y=(gal.Zspec(f));
% plot(x,y,'.','markersize',1); hold on;
[xx,yy,n,s]=densitygrid(x,y,[nbin,nbin],[8,13],[0,0.4]);
[~,h2]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'r');hold on;
% xlim([31,36]);ylim([-3,1]);
xlabel('Mstar');ylabel('z');