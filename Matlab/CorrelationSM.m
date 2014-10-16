%%
xbrd=[5e9,1e10; 1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11; 3.7e11, 6e11];
f=gal.CentralSampleIter>0;
[smmid,zmid2]=skeleton(gal.SMsps(f),gal.Zspec(f),[xbrd(:,1);xbrd(end)]);
f=gal.CentralSampleIter>0&gal.FlagSFR>0&gal.SSFR>10^-1.5;
[smmid,zmid]=skeleton(gal.SMsps(f),gal.Zspec(f),[xbrd(:,1);xbrd(end)]);
%%
% cd /work/Projects/Lensing/outputv4/data/
cd /mnt/charon/Lensing/output/
fmts={'b-','r--','g-','m--'}
files={'WL_GalCenSM1.Active.corr.hdf5','WL_GalCenSM1.Full.corr.hdf5','WL_GalCenSM4.Active.corr.hdf5','WL_GalCenSM4.Full.corr.hdf5'};
masses=[12.01,12.51,12.8,13.11];
mass=masses([1,1,3,3]);
zref=[zmid(1),zmid2(1),zmid(4),zmid2(4)];
%-----------
myfigure;
for i=1:4
file=files{i}
m=10.^mass(i)/1e10;
r=h5read(file,'/shear/seperation');
n=h5read(file,'/shear/numpair');
n=double(n);
nr=h5read(file,'/rand/numpair');
nerr=h5read(file,'/rand/numpair_err');
nmock=h5read(file,'/rand/numrand');
nmock=double(nmock);
rv=1;
% rv=comoving_200b_radius(m,0.3)/1000;
rat=n./(nr.*nmock/max(nmock));
rat_err=nerr./nr;
errorbar(r/rv,rat-1,rat_err,fmts{i});
hold on;
end
% l=legend('L1','L2','L3','L3Tight'); set(l,'interpreter','latex');
plot([1e-2,30],[0,0],'k--')
xlabel('Comoving $R[\rm{Mpc}/h]$');ylabel('$n/n_{rand}\frac{\ }{\ }1$');
set(gca,'xscale','log');
yscale('log');
ylim([1e-2,4]);
xlim([1e-2,3]);
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/GAMASrcCorrelation.eps')
%% WLprof
fmts={'bo','ro','gd','md'};
fmtp={'b','r','g','m'};
files={'WL_GalCenSM1.ActiveBroad.RPhys.Dsig.hdf5','WL_GalCenSM1.PassiveBroad.RPhys.Dsig.hdf5','WL_GalCenSM2.ActiveBroad.RPhys.Dsig.hdf5','WL_GalCenSM2.PassiveBroad.RPhys.Dsig.hdf5'};
masses=[12.01,12.51,12.8,13.11];
%-----------
myfigure;
for i=3:4
file=files{i}
r=h5read(file,'/shear/seperation');
s=h5read(file,'/shear/profile');
cov=h5read(file,'/shear/covariance');
es=sqrt(diag(cov));
s1=h5read(file,'/Star2CenActiveFitted_All/profile');
s2=h5read(file,'/Star2CenFitted_All/profile');
ub=s+es;
lb=s-es; 
% lb(lb<0)=1e-3;
ploterr(r,s,[],{lb,ub},fmts{i},'logx');hold on;
if mod(i,2)
plot(r,s1,'--','color',fmtp{i});
else
plot(r,s2,'--','color',fmtp{i});
end
hold on;
end
% l=legend('L1','L2','L3','L3Tight'); set(l,'interpreter','latex');
% plot([1e-2,30],[0,0],'k--')
xlabel('Physical $R[\rm{Mpc}/h]$');ylabel('Physical $\Delta\Sigma$');
set(gca,'xscale','log');
% yscale('log');
ylim([-2e4,2e4]);
xlim([1e-2,1.3])
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/WLProfHighSMPhysical.eps')