cd /work/Projects/Lensing/data
% Ac=-1.2;An=20.7;Az=2.3;
% Bc=0.94;Bn=-0.67;Bz=0.16;
Ac=2.0;An=17.9;Az=1.5; %r<19.8
Bc=0.65;Bn=-0.50;Bz=0.22;%r<19.8
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
for i=1:numel(gal.CATAID)
    sm=galv6.SM(galv6.CATAID==gal.CATAID(i));
    if isempty(sm)
        sm=NaN;
    end
    gal.SM(i)=sm;
end
clear galv6 DM zgal dm kcorr sm
%% now match grp centrals, replace them with galref
for i=1:numel(grp.Gnum)
    grp.IterCenRef(i)=find(grp.IterCenCATAID(i)==gal.CATAID);
    grp.BCGRef(i)=find(grp.BCGCATAID(i)==gal.CATAID);
end
grp.IterCenSM=gal.SM(grp.IterCenRef);
grp.BCGCenSM=gal.SM(grp.BCGRef);
%% galaxy additional properties
gal.IsIterCen=zeros(size(gal.CATAID));
gal.IsIterCen(grp.IterCenRef)=1;
gal.IsBCG=zeros(size(gal.CATAID));
gal.IsBCG(grp.BCGref)=1;
gal.EnvLevel=(gal.GroupID>0)+gal.IsBCG; %environment: 0: field; 1: sat; 2: bcg.

ind=zeros(max(grp.GroupID),1);
ind(grp.GroupID)=1:numel(grp.GroupID);
gal.LumMass=zeros(size(gal.CATAID));
gal.LumMass(gal.GroupID>0)=grp.LumMass(ind(gal.GroupID(gal.GroupID>0)));
%%
save G3Cv6.mat grp gal
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
%% load the new mocks
for i=1:9
    mock{i}=load(['201310/mockcat_',num2str(i),'.mat']);
end
%% old mocks
% grpmock=fits_load_bintable('201309/DMUG3Cv06/mocks/G3CMockFoFGroupv06.fits',0,1,1);
% galmock=fits_load_bintable('201309/DMUG3Cv06/mocks/G3CMockGalv06.fits',0,1,1);
% grpmock=mvfield(grpmock,{'Nfof','MassProxy','TotFluxProxy'},{'Mult','MassProxyRaw','TotFluxProxyRaw'});
% %Now the calibrated dynmass and luminosity, luminosity mass are called
% %DynMass, Luminosity, LumMass
% grpmock.DynMass=grpmock.MassProxyRaw.*(Ac+An./sqrt(grpmock.Mult)+Az./sqrt(grpmock.Zfof));
% grpmock.Luminosity=grpmock.TotFluxProxyRaw.*(Bc+Bn./sqrt(grpmock.Mult)+Bz./sqrt(grpmock.Zfof));
% [grpmock.LumMass,C3,A3]=luminosity_mass(grpmock);
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