cd /work/Projects/Lensing/data
% Ac=-1.2;An=20.7;Az=2.3;
% Bc=0.94;Bn=-0.67;Bz=0.16;
Ac=2.0;An=17.9;Az=1.5; %r<19.8
% Bc=0.65;Bn=-0.50;Bz=0.22;%r<19.8
Bc=0.086;Bn=-0.37;Bz=0.35; %new calibration, r<19.8
%%
% grp=fits_load_bintable('201308/G3Cv6/G3CFoFGroupv06.fits',0,1,1);
% gal=fits_load_bintable('201308/G3Cv6/G3CGalv06.fits',0,1,1);
%%
grp=fits_load_bintable('201309/DMUG3Cv06/groups/G3CFoFGroupv06.fits',0,1,1);
gal=fits_load_bintable('201309/DMUG3Cv06/groups/G3CGalv06.fits',0,1,1);
grp=mvfield(grp,{'Nfof','MassProxy','TotFluxProxy'},{'Mult','MassProxyRaw','TotFluxProxyRaw'});
%Now the calibrated dynmass and luminosity, luminosity mass are called
%DynMass, Luminosity, LumMass
grp.DynMass=grp.MassProxyRaw.*(Ac+An./sqrt(grp.Mult)+Az./sqrt(grp.Zfof));
grp.Luminosity=grp.TotFluxProxyRaw.*(Bc+Bn./sqrt(grp.Mult)+Bz./sqrt(grp.Zfof));
[grp.LumMass,C3,A3]=luminosity_mass(grp);
%% prepare stellar mass
galv6=fits_load_bintable('/work/Projects/Lensing/data/201308/TilingCatv40_kcorr_z00v03.fits');
dm=@(z) 5*log10(lum_dist(0.3,z)*1e6/10); %DM+5log10(h)
kcorr=@(p1,p2,p3,p4,p5,z) p1.*z.^4+p2.*z.^3+p3.*z.^2+p4.*z+p5;
zgal=galv6.Z_TONRY;
DM=dm(zgal);DM(zgal<=0)=0;
g=galv6.G_MODEL-DM-kcorr(galv6.PCOEFF_G_1,galv6.PCOEFF_G_2,galv6.PCOEFF_G_3,galv6.PCOEFF_G_4,galv6.PCOEFF_G_5,zgal);
i=galv6.I_MODEL-DM-kcorr(galv6.PCOEFF_I_1,galv6.PCOEFF_I_2,galv6.PCOEFF_I_3,galv6.PCOEFF_I_4,galv6.PCOEFF_I_5,zgal);
galv6.SM=10.^(1.15+0.70*(g-i)-0.4*i);%Msun/h^2, since our mag (since DM) is actually Mag-5log10(h)
clear g i
% gal.SM=zeros(size(gal.CATAID));
% for i=1:numel(gal.CATAID)
%     sm=galv6.SM(galv6.CATAID==gal.CATAID(i));
%     if isempty(sm)
%         sm=NaN;
%     end
%     gal.SM(i)=sm;
% end
% clear galv6 DM zgal dm kcorr sm
%% now match grp centrals, replace them with galref
grp.IterCenRef=zeros(size(grp.GroupID));
grp.BCGRef=zeros(size(grp.GroupID));
grp.IterCenSM=zeros(size(grp.GroupID));
grp.BCGSM=zeros(size(grp.GroupID));
for i=1:numel(grp.GroupID)
    grp.IterCenRef(i)=find(grp.IterCenCATAID(i)==gal.CATAID);
    grp.BCGRef(i)=find(grp.BCGCATAID(i)==gal.CATAID);
end
grp.IterCenSM=gal.SM(grp.IterCenRef);
grp.BCGSM=gal.SM(grp.BCGRef);
%% galaxy additional properties
gal.IsIterCen=zeros(size(gal.CATAID));
gal.IsIterCen(grp.IterCenRef)=1;
gal.IsBCG=zeros(size(gal.CATAID));
gal.IsBCG(grp.BCGRef)=1;
gal.EnvLevel=(gal.GroupID>0)+gal.IsBCG; %environment: 0: field; 1: sat; 2: bcg.

ind=zeros(max(grp.GroupID),1);
ind(grp.GroupID)=1:numel(grp.GroupID);
gal.LumMass=zeros(size(gal.CATAID));
gal.LumMass(gal.GroupID>0)=grp.LumMass(ind(gal.GroupID(gal.GroupID>0)));
%%
save G3Cv6.mat grp gal
%% prepare galaxy_number_density-redshift relation
data=importdata('LumFunc/lf_gama.cat.2.0.01_z_0.51_swml_z.0.05.0.10.best.txt',' ',2);
lumfunc_data=data.data;
dm=@(z) 5*log10(lum_dist(0.3,z)*1e6/10); %DM+5log10(h)
pkcorr=[2.3843,3.5902,0.5237,1.0226,0.2085]; %k correction polynomial
kecorr=@(z) polyval(pkcorr,z-0.2)+1.75*z;
Mlim=@(z) 19.8-dm(z)-kecorr(z);
cols=5:9;
z=0.01:0.02:0.52;
zbin=[0:0.1:0.49,0.6];
x=lumfunc_data(:,1);
x(x==0)=lumfunc_data(x==0,2);
[~,ibins]=histc(z,zbin);
Nav=z;
for i=1:numel(z)
%     y=lumfunc_data(:,cols(ibins(i)));
    y=lumfunc_data(:,3);
    fun = fit(x,y,'linearinterp');
    Nav(i)=integrate(fun,Mlim(z(i)),-36);
end
figure;
loglog(z,Nav,'.');hold on;
plot(z-0.005,pchip(z,Nav,z-0.005),'r-');
Nav=[z',Nav'];
% save LumFunc/GAMA-Ndens-z-global.txt Nav -ascii
%%
fun = fit(x,lumfunc_data(:,3),'linearinterp');
Nav_abs=integrate(fun,-13,-26); %=0.16, average number density within M=[-26,-13]mag, multiply this to convert VolMult to Nabs.
%%
h=0.7;
dm=@(z) 5*log10(lum_dist(0.3,z)*1e6/h/10);
sfz=@(x,a) -(exp(-x)./(a+1).*x.^(a+1))+gammainc(x,2+a,'upper')*gamma(2+a)/(a+1);
ffz2=@(z,ml,Mc,a) sfz(10.^(-0.4*(ml-Mc-dm(z))),a)./sfz(10.^(-0.4*(-13-Mc)),a); %corrected to M=-13, slightly fainter than the limiting m=19.4 at z=0.01 is M=-13.78
Mc=-20.73+5*log10(h);
a=-1.26;
loveday_nrat=@(z) 1./ffz2(z,19.8,Mc-0.7*z,a);
Nav=load('LumFunc/GAMA-Ndens-z.txt','-ascii');
Nav2=load('LumFunc/GAMA-Ndens-z-global.txt','-ascii');
figure;
semilogy(Nav(:,1),0.16./Nav(:,2),'r');
hold on;
semilogx(Nav2(:,1),0.16./Nav2(:,2),'g');
semilogx(Nav(:,1),loveday_nrat(Nav(:,1)),'b');
%%
Nav=load('LumFunc/GAMA-Ndens-z.txt','-ascii');
Nav2=load('LumFunc/GAMA-Ndens-z-global.txt','-ascii');
figure;
semilogy(Nav(:,1),Nav(:,2),'r');
hold on;
plot(Nav2(:,1),Nav2(:,2),'g');
figure;
semilogx(Nav(:,1),Nav(:,2)./Nav2(:,2),'r');
%% extract number-density directly
declim=[-2,3;-3,2;-2,3];
Area=12*pi/180*sum(sind(declim(:,2))-sind(declim(:,1)));
zbin=[0:0.03:0.52]';
% zbin=logspace(-2,log10(0.5),10)';
zmid=(zbin(1:end-1)+zbin(2:end))/2;
rbin=comoving_dist(0.3,0.7,zbin);
N=histc(gal.Z,zbin);
n=N(1:end-1)./(Area/3*diff(rbin.^3)); %comoving density
en=sqrt(N(1:end-1))./(Area/3*diff(rbin.^3));
Ndens=[zmid,n,en];
save GAMAII-Ndens-z.txt Ndens -ascii
%% compare the number densities
figure;
semilogy(zmid,n,'rx')
hold on;
plot(zmid,0.1*exp(-10.4*zmid),'r-');
z=0:0.01:0.6;
plot(z,pchip(zmid,n,z),'r--');
plot(Nav2(:,1),Nav2(:,2),'g');
semilogx(Nav(:,1),0.16./loveday_nrat(Nav(:,1)),'b');
legend('Data','0.1*exp(-10.4*z)','PederGlobal','Loveday');
%% add multiplicity volume=N/n, comoving volume in Mpc/h
load G3Cv6
% Nav=load('LumFunc/GAMA-Ndens-z-global.txt','-ascii');
% Nav=load('GAMAII-Ndens-z.txt','-ascii');
% avNdens=@(z) pchip(Nav(:,1),Nav(:,2),z);
% grp.VolMult=grp.Mult./avNdens(grp.Zfof);
grp.Luminosity=grp.TotFluxProxyRaw.*(Bc+Bn./sqrt(grp.Mult)+Bz./sqrt(grp.Zfof));
[grp.LumMass,C3,A3]=luminosity_mass(grp);
save G3Cv6 grp -append
%%
% myfigure;
% semilogx(grp.Mult,grp.MassAfunc./grp.DynMass,'.');
% xlabel('Multiplicity');ylabel('MassAfunc/(MassProxy*Afunc)');
% print('-depsc','/work/Projects/Lensing/outputv4/MassAfuncCheck.eps');
% %%
% myfigure;
% semilogx(grp.Mult,grp.LumMassBfunc./grp.Luminosity,'.');
% xlabel('Multiplicity');ylabel('LumBfunc/(FluxProxy*Bfunc)');
% print('-depsc','/work/Projects/Lensing/outputv4/LumBfuncCheck.eps');
% figure;
% semilogx(grp.TotFluxProxyRaw,grp.LumMassB./grp.TotFluxProxyRaw,'.');
%%
% grp.DynMass=grp.MassAfunc;
% grp.Luminosity=grp.LumMassBfunc; %this column is problematic. do not trust it.
%% old mocks
grpmock=fits_load_bintable('201309/DMUG3Cv06/mocks/G3CMockFoFGroupv06.fits',0,1,1);
galmock=fits_load_bintable('201309/DMUG3Cv06/mocks/G3CMockGalv06.fits',0,1,1);
grpmock=mvfield(grpmock,{'Nfof','MassProxy','TotFluxProxy'},{'Mult','MassProxyRaw','TotFluxProxyRaw'});
%Now the calibrated dynmass and luminosity, luminosity mass are called
%DynMass, Luminosity, LumMass
grpmock.DynMass=grpmock.MassProxyRaw.*(Ac+An./sqrt(grpmock.Mult)+Az./sqrt(grpmock.Zfof));
grpmock.Luminosity=grpmock.TotFluxProxyRaw.*(Bc+Bn./sqrt(grpmock.Mult)+Bz./sqrt(grpmock.Zfof));
[grpmock.LumMass,C3,A3]=luminosity_mass(grpmock);
Nav=load('GAMAII-Ndens-z.txt','-ascii');
avNdens=@(z) pchip(Nav(:,1),Nav(:,2),z);
grpmock.VolMult=grpmock.Mult./avNdens(grpmock.Zfof);
%% split the mock
grpmockall=grpmock;
galmockall=galmock;
for vol=1:9
    names=fieldnames(grpmockall);
    for i=1:numel(names)
        grpmock.(names{i})=grpmockall.(names{i})(grpmockall.Volume==vol);
    end
    names=fieldnames(galmockall);
    for i=1:numel(names)
        galmock.(names{i})=galmockall.(names{i})(galmockall.Volume==vol);
    end
    save(['201309/mockcat_',num2str(vol),'.mat'],'grpmock','galmock');
end
%% now match grp centrals, replace them with galref, and append halomass
for vol=1:9
    load(['201309/mockcat_',num2str(vol),'.mat'])
    grpmock.IterCenRef=zeros(size(grpmock.GroupID));
    grpmock.BCGRef=zeros(size(grpmock.GroupID));
    for i=1:numel(grpmock.GroupID)
        grpmock.IterCenRef(i)=find(grpmock.IterCenGalID(i)==galmock.GalID);
        grpmock.BCGRef(i)=find(grpmock.BCGGalID(i)==galmock.GalID);
    end
    grpmock.MIter=galmock.HaloMass(grpmock.IterCenRef);
    grpmock.MBCG=galmock.HaloMass(grpmock.BCGRef);
    save(['201309/mockcat_',num2str(vol),'.mat'],'grpmock','-append');
end
%% calculate VolMult with mocks' own dn/dv
declim=[-1,3;-2,2;-2,2];Area=12*pi/180*sum(sind(declim(:,2))-sind(declim(:,1)));%note the area is G3Cv4 area!
zbin=[0:0.03:0.52]';
zmid=(zbin(1:end-1)+zbin(2:end))/2;
rbin=comoving_dist(0.3,0.7,zbin);
Vbin=Area/3*diff(rbin.^3);
for i=1:9
    load(['201309/mockcat_',num2str(i),'.mat']);
    Nmock=histc(galmock.Z,zbin);
    nmock=Nmock(1:end-1)./Vbin; %comoving density
    grpmock.VolMult=grpmock.Mult./pchip(zmid,nmock,grpmock.Zfof);% calibrate with its own dN/dV
    save(['201309/mockcat_',num2str(i),'.mat'],'grpmock','-append');
end
%% new lum calibration
Bc=0.086;Bn=-0.37;Bz=0.35; %new calibration, r<19.8
for i=1:9
    load(['201309/mockcat_',num2str(i),'.mat']); 
    grpmock.Luminosity=grpmock.TotFluxProxyRaw.*(Bc+Bn./sqrt(grpmock.Mult)+Bz./sqrt(grpmock.Zfof)); %redo the calibration
    [grpmock.LumMass,C1,A1]=luminosity_mass(grpmock);
    save(['201309/mockcat_',num2str(i),'.mat'],'grpmock','-append');
end
%%
myfigure;
loglog(galmockall{1}.Mstar,galmockall{1}.Mcolor,'.');
hold on;
plot([1e6,1e13],[1e6,1e13],'k-');
xlabel('$Mstar$');ylabel('$Mcolor$');
print('-dpng','/work/Projects/Lensing/outputv4/MockMstar-Mcolor.png');
%%
% figure;
[xm,ym,yl,xmean,ymean,ysig]=skeleton(galmockall{1}.Mcolor(galmockall{1}.is_central>0),log10(galmockall{1}.Mhalo(galmockall{1}.is_central>0)),logspace(10,12,10),0.683);
plot(xm,ym,'y-');
hold on;
plot(xm,yl(:,1),'y--');
plot(xm,yl(:,2),'y--');
f=grpmockall{1}.Mult>5;
[xm,ym,yl,xmean,ymean,ysig]=skeleton(grpmockall{1}.McolorIter(f),log10(grpmockall{1}.MIter(f)),logspace(10,12,10),0.683);
plot(xm,ym,'b-');
hold on;
plot(xm,yl(:,1),'b--');
plot(xm,yl(:,2),'b--');