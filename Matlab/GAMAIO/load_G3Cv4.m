cd /work/Projects/Lensing/data
grpmock=loadGAMAcsv('G3Cv4/mock/G3CMock1FoFGroup194v04.dat',4);
halomock=loadGAMAcsv('G3Cv4/mock/G3CMock1HaloGroup194v04.dat',4);
galmock=loadGAMAcsv('G3Cv4/mock/G3CMock1Galv04.dat',4);
refmock=loadGAMAcsv('G3Cv4/mock/G3CMock1Ref194v04.dat',4);
Ac=-1.2;An=20.7;Az=2.3;
Bc=0.94;Bn=-0.67;Bz=0.16;
halomock.MassProxy=halomock.MassProxy.*(Ac+An./sqrt(halomock.Mult)+Az./sqrt(halomock.Zfof));
% halomock.MassProxy=halomock.MassProxy.*dynmass_calibration(halomock); %hardly any difference, when using the tabulated corrections
halomock.TotFluxProxy=halomock.TotFluxProxy.*(Bc+Bn./sqrt(halomock.Mult)+Bz./sqrt(halomock.Zfof));
[halomock.LumMass,C1,A1]=luminosity_mass(halomock);
%[grpmock.DynMass,CD,AD]=dynamical_mass(grpmock,mockmass.MIter);
grpmock.MassProxy=grpmock.MassProxy.*(Ac+An./sqrt(grpmock.Mult)+Az./sqrt(grpmock.Zfof));
% grpmock.MassProxy=grpmock.MassProxy.*dynmass_calibration(grpmock);
tmp=grpmock.TotFluxProxy;
grpmock.TotFluxProxy=grpmock.TotFluxProxy.*(Bc+Bn./sqrt(grpmock.Mult)+Bz./sqrt(grpmock.Zfof));
grpmock.TotFluxProxyRaw=tmp;
clear tmp;
[grpmock.LumMass,C2,A2]=luminosity_mass(grpmock);

grp=loadGAMAcsv('G3Cv4/group/G3CFoFGroup194v04.dat',4);
ref=loadGAMAcsv('G3Cv4/group/G3CRef194v04.dat',4);
grp.MassProxyRaw=grp.MassProxy;
grp.TotFluxProxyRaw=grp.TotFluxProxy;
grp.MassProxy=grp.MassProxy.*(Ac+An./sqrt(grp.Mult)+Az./sqrt(grp.Zfof));
grp.TotFluxProxy=grp.TotFluxProxy.*(Bc+Bn./sqrt(grp.Mult)+Bz./sqrt(grp.Zfof));
[grp.LumMass,C3,A3]=luminosity_mass(grp);
grp.HybMass=sqrt(grp.MassProxy.*grp.LumMass);

gal=loadGAMAcsv('G3Cv4/group/G3CGalv04.dat',4);

linkmock=[refmock.GalID,refmock.GroupID];  %galID: index (row number) for galaxies in the galaxy file
linkmock=sortrows(linkmock,2);
link=[ref.GalID,ref.GroupID];
link=sortrows(link,2);

%% mockgal additional properties
dir='/work/Projects/Lensing/data/2011_02';
Tz0=importdata([dir,'/Tz0.new'],' ',1);
galmock.HaloMass=Tz0.data(:,8);
galmock.Mag=Tz0.data(:,5);
%%
GrpOffset=sum(~linkmock(:,2)); %skip isolated galaxies
a=[];
for gid=1:10
    x=galmock.RA(linkmock(GrpOffset+(1:grpmock.Mult(gid)),1));
    y=galmock.DEC(linkmock(GrpOffset+(1:grpmock.Mult(gid)),1));
    r=acos(ccdist([x,y],[grpmock.IterCenRA(gid),grpmock.IterCenDEC(gid)]));
    r100=max(r)*AD_dist_flat(0.3,0,grpmock.Zfof(gid)).*(1+grpmock.Zfof(gid));%comoving
    a=[a;r100,grpmock.Rad100(gid)];%the Rad*** are comoving!
    GrpOffset=GrpOffset+grpmock.Mult(gid);
end
%%  mark-out un-matched halos?
% id2ind=zeros(max(max(halomock.HaloID),max(galmock.HaloID)),1);
% id2ind(halomock.HaloID)=1:numel(halomock.HaloID);
% galmock.HaloMass(~id2ind(galmock.HaloID))=1e9; %use 1e9 as NaN here.
% clear id2ind
%% galmock halo mass
mockmass.Complexity=grpmock.GroupID;  %number of different group masses
mockmass.Mmean=grpmock.GroupID;   % average mass
mockmass.Gmean=grpmock.GroupID;   % geometric average mass, or average in logspace
mockmass.Mmode=grpmock.GroupID;   % most frequent mass
mockmass.Purity=grpmock.GroupID;  % frequency of Mmode
mockmass.MBCG=galmock.HaloMass(grpmock.BCGRef); %group mass for BCG
mockmass.MIter=galmock.HaloMass(grpmock.IterCenRef); %group mass for BCG
GrpOffset=sum(~linkmock(:,2)); %skip isolated galaxies
for gid=1:numel(grpmock.GroupID)
    mass=galmock.HaloMass(linkmock(GrpOffset+(1:grpmock.Mult(gid)),1));
    mag=galmock.Mag(linkmock(GrpOffset+(1:grpmock.Mult(gid)),1));
    ids=galmock.HaloID(linkmock(GrpOffset+(1:grpmock.Mult(gid)),1));
    ids=[ids,galmock.GalID(linkmock(GrpOffset+(1:grpmock.Mult(gid)),1))];
    GrpOffset=GrpOffset+grpmock.Mult(gid);
    ids=sortrows(ids);
    mockmass.Complexity(gid)=sum(logical(diff(ids(:,1))))+1;
    mockmass.Mmean(gid)=mean(mass);
    mockmass.Gmean(gid)=geomean(mass(logical(mass)));  %geomean of non-zero elements
    [mockmass.Mmode(gid),mockmass.Purity(gid)]=mode(ids(:,1));
    mid=ids(ids(:,1)==mockmass.Mmode(gid),2);mid=mid(1); %take GalID of the mode galaxy
    mockmass.Mmode(gid)=galmock.HaloMass(mid);
    if mockmass.Purity(gid)==1
%         mockmass.Mmode(gid)=max(mass);
        [~,minmagid]=min(mag);
        mockmass.Mmode(gid)=mass(minmagid);
    end
end
mockmass.Purity=mockmass.Purity./grpmock.Mult;
clear Tz0 dir GrpOffset mass mag gid ids mid minmagid
grpmock.MIter=mockmass.MIter;
%% gal additional properties from G3Cv2
g09=loadGAMAcsv('G3Cv2/group/cutcatG09.csv');
g12=loadGAMAcsv('G3Cv2/group/cutcatG12.csv');
g15=loadGAMAcsv('G3Cv2/group/cutcatG15.csv');
galv2=structcat([g09,g12,g15]);
clear g09 g12 g15
gal.SM=galv2.STELLAR_MASS;
gal.Upetro_SDSS=galv2.PETROMAG_U_SDSS;
gal.Rpetro_SDSS=galv2.PETROMAG_R_SDSS;
gal.A_u=galv2.A_u;
gal.A_r=galv2.A_r;
gal.AB_u=galv2.AB_u;
gal.AB_r=galv2.AB_r;
gal.SM(gal.SM==9999)=NaN;
clear galv2
%%
galv6=fits_load_bintable('/work/Projects/Lensing/data/201308/TilingCatv40_kcorr_z00v03.fits');
hubble=0.7;
dm=@(z) 5*log10(lum_dist(0.3,z)*1e6/hubble/10);
kcorr=@(p1,p2,p3,p4,p5,z) p1.*z.^4+p2.*z.^3+p3.*z.^2+p4.*z+p5;
zgal=galv6.Z_TONRY;
DM=dm(zgal);DM(zgal<=0)=0;
g=galv6.G_MODEL-DM-kcorr(galv6.PCOEFF_G_1,galv6.PCOEFF_G_2,galv6.PCOEFF_G_3,galv6.PCOEFF_G_4,galv6.PCOEFF_G_5,zgal);
i=galv6.I_MODEL-DM-kcorr(galv6.PCOEFF_I_1,galv6.PCOEFF_I_2,galv6.PCOEFF_I_3,galv6.PCOEFF_I_4,galv6.PCOEFF_I_5,zgal);
galv6.SM=10.^(1.15+0.70*(g-i)-0.4*i)*hubble^2;%Msun/h^2
clear g i
% gal.SMnew=gal.SM;
gal.SMsps=gal.SM;
for i=1:numel(gal.A_u)
    sm=galv6.SM(galv6.CATAID==gal.GalID(i));
    if isempty(sm)
        sm=NaN;
    end
    gal.SM(i)=sm;
    sm=galv6.logSM_sps(galv6.CATAID==gal.GalID(i)&galv6.flag_sps>0&galv6.Z_TONRY>0.003&galv6.NQ>2);
    if isempty(sm)
        sm=NaN;
    end
    gal.SMsps(i)=10.^sm;
end
clear galv6 DM zgal dm kcorr sm
%% SFR
dir='/work/Projects/Lensing/data/additional/';
% SFRtable=fits_load_bintable([dir,'catemis_v1.fits']);
% SFRtable=mvfield(SFRtable,'CATA_INDEX_1','GalID');
% id2ind=zeros(max(max(gal.GalID),max(SFRtable.GalID)),1);
% id2ind(SFRtable.GalID)=1:numel(SFRtable.GalID);
% gal.logSFR=zeros(size(gal.GalID))+NaN;
% gal.logSSFR=gal.logSFR;
% ind=id2ind(gal.GalID);
% f=logical(ind);
% gal.logSFR(f)=SFRtable.Log_SFR(ind(f));
% gal.logSSFR(f)=SFRtable.Log_SSFR(ind(f));
% dlmwrite([dir,'logSFR.dat'],[gal.logSFR,gal.logSSFR]);
% clear id2ind ind f dir
SFR=dlmread([dir,'logSFR.dat']);
gal.logSFR=SFR(:,1);
gal.logSSFR=SFR(:,2);
clear dir SFR
% mid value: 0.7
%%
hosts=zeros(max(gal.GalID),1);
hosts(ref.GalID)=ref.GroupID;

cenhosts=zeros(max(gal.GalID),1);
cenhosts(grp.IterCenRef)=1;
cenhosts(logical(cenhosts))=hosts(logical(cenhosts));

bcghosts=zeros(max(gal.GalID),1);
bcghosts(grp.BCGRef)=1;

f=zeros(numel(gal.GalID),1);
f(logical(hosts(gal.GalID)))=1;

f(logical(bcghosts(gal.GalID)))=2;
gal.IsIterCen=cenhosts(gal.GalID);
gal.EnvLevel=f; %environment: 0: field; 1: sat; 2: bcg.

gal.GroupID=hosts(gal.GalID);
clear f hosts cenhosts bcghosts

ind=zeros(max(grp.GroupID),1);
ind(grp.GroupID)=1:numel(grp.GroupID);
gal.LumMass=zeros(size(gal.GalID));
gal.LumMass(gal.GroupID>0)=grp.LumMass(ind(gal.GroupID(gal.GroupID>0)));
%%
SM=zeros(max(gal.GalID),1);
SM(gal.GalID)=gal.SM; %Msun/h^2
grp.IterCenSM=SM(grp.IterCenRef);
grp.BCGSM=SM(grp.BCGRef);
GrpSM=zeros(max(grp.GroupID),1);
for i=1:numel(ref.GalID)
    if ref.GroupID(i)>0
    GrpSM(ref.GroupID(i))=GrpSM(ref.GroupID(i))+SM(ref.GalID(i));
    end
end
grp.TotSM=GrpSM(grp.GroupID);
%% add stellar-population-synthesis mass
SM=zeros(max(gal.GalID),1);
SM(gal.GalID)=gal.SMsps; %Msun/h^2
grp.IterCenSMsps=SM(grp.IterCenRef);
grp.BCGSMsps=SM(grp.BCGRef);
%%

MvInMvb=0.893;
MvcInMvb=0.733;
halo2star=@(M,A,M0,alpha,beta,gamma) A./((M/M0).^alpha+(M/M0).^beta).^gamma.*M;
halo2starP=@(M,P) halo2star(M,P(1),P(2),P(3),P(4),P(5)); %Msun/h^2 for Mstar, Msun/h for Mhalo
% Pwang=[2*10^10.17/3.21e11*0.73^2*MvcInMvb,3.21e11/MvcInMvb,1-2.42,1-0.29,1];%unified, close to ling
Pwang=[2*10^10.23/3.43e11*0.73^2*MvcInMvb,3.43e11/MvcInMvb,1-2.56,1-0.34,1];%DR7, 2013
Pguo=[0.129*0.73*MvcInMvb,10^11.4*0.73/MvcInMvb,-0.926,0.261,2.440];  % M200c, convert to M200b
Pmoster=[2*0.0282*0.72*MvInMvb,10^11.884*0.72/MvInMvb,-1.057,0.556,1]; % tophat Mvir, convert to M200b
Pling=[2*10^-1.73,10^11.70,-1.16,0.71,1];%M200b, Wang Lingyu 2013.

M=logspace(9,17,100);

Msyang=10^10.86*(M/10^12.08).^(0.22+1.61)./(1+M/10^12.08).^1.61; %M180b, not yet converted
Mswang=halo2starP(M,Pwang);
Msguo=halo2starP(M,Pguo);
Msmoster=halo2starP(M,Pmoster);
Msling=halo2starP(M,Pling);

grp.HODMassIterMoster=spline(Msmoster,M,grp.IterCenSM);
% grp.HODMassBCG=spline(Msmoster,M,grp.BCGSM);
grp.HODMassIterGuo=spline(Msguo,M,grp.IterCenSM);
grp.HODMassIterYang=spline(Msyang,M,grp.IterCenSM);
grp.HODMassIterLing=spline(Msling,M,grp.IterCenSM);
grp.HODMassIterWang=spline(Mswang,M,grp.IterCenSM);
clear MvInMvb MvcInMvb Pwang Pguo Pmoster Pling M Mswang Msguo Msmoster Msling Msyang
%%
% save GAMAgalv4up6 gal -v7.3 %updated with G3Cv6 StellarMass and HODmass
% save GAMAgrpv4up6 grp -v7.3
%%

cd /work/Projects/Lensing/outputv4