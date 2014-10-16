clear
clc
cd /work/Projects/Lensing/data
% Ac=-1.2;An=20.7;Az=2.3;
% Bc=0.94;Bn=-0.67;Bz=0.16;
% Ac=2.0;An=17.9;Az=1.5; %r<19.8
% Bc=0.65;Bn=-0.50;Bz=0.22;%r<19.8
% Bc=0.086;Bn=-0.37;Bz=0.35; %new calibration, r<19.8
%%
load GAMAgrpv4up6
load GAMAgalv4up6
%%
grp=mvfield(grp,{'MassProxy','TotFluxProxy'},{'DynMass','Luminosity'});
%Now the calibrated dynmass and luminosity, luminosity mass are called
%DynMass, Luminosity, LumMass
%% add new luminosity and lummass calibration (trial, may not be useful at all)
Bc=0.086;Bn=-0.37;Bz=0.35; %new calibration, r<19.8
tmp=grp.Luminosity;
grp.Luminosity=grp.TotFluxProxyRaw.*(Bc+Bn./sqrt(grp.Mult)+Bz./sqrt(grp.Zfof));
[grp.LumMass2,C3,A3]=luminosity_mass(grp); 
grp.Luminosity2=grp.Luminosity; % rename
grp.Luminosity=tmp; % restore
clear tmp
%% extract number-density directly
declim=[-1,3;-2,2;-2,2];
Area=12*pi/180*sum(sind(declim(:,2))-sind(declim(:,1)));
zbin=[0:0.03:0.52]';
zmid=(zbin(1:end-1)+zbin(2:end))/2;
rbin=comoving_dist(0.3,0.7,zbin);
Vbin=Area/3*diff(rbin.^3);
N=histc(gal.Zspec,zbin);
n=N(1:end-1)./Vbin; %comoving density, (Mpc/h)^-3
en=sqrt(N(1:end-1))./Vbin;
Ndens=[zmid,n,en];
save GAMAI-Ndens-z.txt Ndens -ascii
%% VolMult
% Ndens=load('GAMAI-Ndens-z.txt','-ascii');
avNdens=@(z) pchip(Ndens(:,1),Ndens(:,2),z);
grp.VolMult=grp.Mult./avNdens(grp.Zfof);
%%
save G3Cv4up8/G3Cv4up8.mat grp gal
%% mocks: Append VolMult, and restore 19.4 lum calibration
Ac=-1.2;An=20.7;Az=2.3;
Bc=0.94;Bn=-0.67;Bz=0.16;
% Nav=load('GAMAI-Ndens-z.txt','-ascii');
% avNdens=@(z) pchip(Nav(:,1),Nav(:,2),z);

declim=[-11.25,-2.75;-2,4;-3,3;-3,3];dra=[11.5;14;14;14];Area=pi/180*sum((sind(declim(:,2))-sind(declim(:,1))).*dra);%much larger area!!
zbin=[0:0.03:0.52]';
zmid=(zbin(1:end-1)+zbin(2:end))/2;
rbin=comoving_dist(0.3,0.7,zbin);
Vbin=Area/3*diff(rbin.^3);

for i=1:9
    load(['201310/mockcat_',num2str(i),'.mat']);
    grpmock.DynMass=grpmock.MassProxyRaw.*(Ac+An./sqrt(grpmock.Mult)+Az./sqrt(grpmock.Zfof));
    [grpmock.LumMass2,C3,A3]=luminosity_mass(grpmock); %r<19.8's new calibration, recalibrate to old dynmass for consistency
    grpmock.Luminosity2=grpmock.Luminosity; % r<19.8's new calibration
    grpmock.Luminosity=grpmock.TotFluxProxyRaw.*(Bc+Bn./sqrt(grpmock.Mult)+Bz./sqrt(grpmock.Zfof)); %standard calibration
    [grpmock.LumMass,C3,A3]=luminosity_mass(grpmock); %standard calibration
    % trim galaxies
    names=fieldnames(galmock);
    f=galmock.r_mag<=19.4;
    for j=1:numel(names)
        galmock.(names{j})=galmock.(names{j})(f);
    end
    % Ndens and Vmult
    Nmock=histc(galmock.Z,zbin);
    nmock=Nmock(1:end-1)./Vbin; %comoving density
    figure;semilogy(zmid,n,'r',zmid,nmock,'g');
    grpmock.VolMult=grpmock.Mult./pchip(zmid,nmock,grpmock.Zfof);% calibrate with its own dN/dV
%   grpmock.VolMult=grpmock.Mult./avNdens(grpmock.Zfof);
    
    save(['G3Cv4up8/mockcat_',num2str(i),'.mat'],'grpmock','galmock');
end
%% (save for Lingyu: Caution: the galaxy file has been trimmed and so does the galaxy indices!)
for i=1:9
    load(['G3Cv4up8/mockcat_',num2str(i),'.mat']);
    fp=fopen(['G3Cv4up8/Lingyu/grpmore_',num2str(i),'.dat'],'w');
    fprintf(fp,'Gnum,IterIndex,BCGIndex,DynMass_Msunh,Luminosity_Lsunh2,LumMass_Msunh\n');
    for j=1:numel(grpmock.Gnum)
        fprintf(fp,'%d,%d,%d,%g,%g,%g\n',grpmock.Gnum(j),grpmock.IterRef(j),grpmock.BCGRef(j),grpmock.DynMass(j),grpmock.Luminosity(j),grpmock.LumMass(j));
    end
    fclose(fp);
end
%% systematic-corrected WL mass
HybMassZ=@(sigma,L,z) 10.^(11.15-0.14)*sigma.^(1.31+0.03).*(L/1e12).^(0.78+0.02).*(1+z).^(-5.79+0.18);
DynVol=@(sigma,V) 10.^(8.49-0.13)*V.^(0.61+0.02).*sigma.^1.28;
LumNew=@(L) 10.^(14.23-0.07)*(L/1e12).^(1.08+0.01);
LumVol=@(L,V) 10.^(17.46-0.40)*(L/1e12).^(1.99-0.10).*V.^(-0.92+0.10);
grp.WLMassVdLumZ=HybMassZ(grp.VelDisp,grp.Luminosity,grp.Zfof);
grp.WLMassVdVn=DynVol(grp.VelDisp,grp.VolMult);
grp.WLMassLum=LumNew(grp.Luminosity);
grp.WLMassLumVn=LumVol(grp.Luminosity,grp.VolMult);
save('G3Cv4up8/G3Cv4up8.mat','grp','-append');
%% save real group for Lingyu
load G3Cv4up8/G3Cv4up8.mat
fp=fopen(['G3Cv4up8/Lingyu/grpmore_GAMA-I.dat'],'w');
fprintf(fp,'Gnum,IterIndex,BCGIndex,DynMass/(Msun/h),Luminosity/(Lsun/h^2),LumMass/(Msun/h)\n');
for j=1:numel(grp.GroupID)
    fprintf(fp,'%d,%d,%d,%g,%g,%g\n',grp.GroupID(j),grp.IterCenRef(j),grp.BCGRef(j),grp.DynMass(j),grp.Luminosity(j),grp.LumMass(j));
end
fclose(fp);
%% save WL mass
load G3Cv4up8/G3Cv4up8.mat
fp=fopen(['G3Cv4up8/Lingyu/grpmore_GAMA-I_Lensing.dat'],'w');
fprintf(fp,'VolMult/(Mpc/h)^3,WLMassVdLumZ/(Msun/h),WLMassLum/(Msun/h),WLMassLumVn/(Msun/h),WLMassVdVn/(Msun/h)\n');
for j=1:numel(grp.GroupID)
    fprintf(fp,'%g,%g,%g,%g,%g\n',grp.VolMult(j),grp.WLMassVdLumZ(j),grp.WLMassLum(j),grp.WLMassLumVn(j),grp.WLMassVdVn(j));
end
fclose(fp);
%% save mock WL mass
for i=1:9
    load(['G3Cv4up8/mockcat_',num2str(i),'.mat']);
    grpmock.WLMassVdLumZ=HybMassZ(grpmock.VelDisp,grpmock.Luminosity,grpmock.Zfof);
    grpmock.WLMassVdVn=DynVol(grpmock.VelDisp,grpmock.VolMult);
    grpmock.WLMassLum=LumNew(grpmock.Luminosity);
    grpmock.WLMassLumVn=LumVol(grpmock.Luminosity,grpmock.VolMult);
    fp=fopen(['G3Cv4up8/Lingyu/grpmore_WLmass_',num2str(i),'.dat'],'w');
    fprintf(fp,'VolMult/(Mpc/h)^3,WLMassVdLumZ/(Msun/h),WLMassLum/(Msun/h),WLMassLumVn/(Msun/h),WLMassVdVn/(Msun/h)\n');
    for j=1:numel(grpmock.Gnum)
        fprintf(fp,'%g,%g,%g,%g,%g\n',grpmock.VolMult(j),grpmock.WLMassVdLumZ(j),grpmock.WLMassLum(j),grpmock.WLMassLumVn(j),grpmock.WLMassVdVn(j));
    end
    fclose(fp);
end
%%
for i=1:9
    load(['G3Cv4up8/mockcat_',num2str(i),'.mat']);
    grpmock.WLMassVdLumZ=HybMassZ(grpmock.VelDisp,grpmock.Luminosity,grpmock.Zfof);
    grpmock.WLMassVdVn=DynVol(grpmock.VelDisp,grpmock.VolMult);
    grpmock.WLMassLum=LumNew(grpmock.Luminosity);
    grpmock.WLMassLumVn=LumVol(grpmock.Luminosity,grpmock.VolMult);
    save(['G3Cv4up8/mockcat_',num2str(i),'.mat'],'grpmock','-append');
end
%%
load G3Cv6.mat
fp=fopen(['G3Cv4up8/Lingyu/grpmore_GAMA-II.dat'],'w');
fprintf(fp,'Gnum,IterIndex,BCGIndex,DynMass_Msunh,Luminosity_Lsunh2,LumMass_Msunh\n');
for j=1:numel(grp.GroupID)
    fprintf(fp,'%d,%d,%d,%g,%g,%g\n',grp.GroupID(j),grp.IterCenRef(j),grp.BCGRef(j),grp.DynMass(j),grp.Luminosity(j),grp.LumMass(j));
end
fclose(fp);
%% add SFR, and correct index ref
for i=1:9
    mock=load(['201310/mockcat_',num2str(i),'.mat']);
    load(['G3Cv4up8/mockcat_',num2str(i),'.mat']);
%     galmock.SFR=mock.galmock.SFR(mock.galmock.r_mag<=19.4);
    newInd=zeros(size(mock.galmock.g));
    newInd(mock.galmock.r_mag<=19.4)=1:numel(galmock.g);
    grpmock.IterRef=newInd(grpmock.IterRef);
    if any(grpmock.IterRef==0)
        error('IterRef error');
    end
%     grpmock.SFRIter=mock.grpmock.SFRIter;
%     if numel(galmock.SFR)~=numel(galmock.g)
%         error('SFR array shape not correct');
%     end
    save(['G3Cv4up8/mockcat_',num2str(i),'.mat'],'grpmock','galmock');
end
%% 2df selection
clear
clc
ngp=importdata('2dF/redcat.cal0b.genuine.2df_ngp',' ',2);
sgp=importdata('2dF/redcat.cal0b.genuine.2df_sgp',' ',2);
twodf=[ngp.data;sgp.data];

n=twodf(:,8)./twodf(:,5);
z=twodf(:,3);
figure; semilogy(z,n,'.')
hold on;
[xm,ym,yl,xmean,ymean,yerr]=skeleton(z,n,0:0.01:0.3,1);
plot(xm,ym,'r');
plot(xm,yl,'r--');
plot(xm,ymean,'ko');
Ndens=load('GAMAI-Ndens-z.txt','-ascii');
plot(Ndens(:,1),Ndens(:,2),'g');
twodfnz=@(z) pchip(xm,ym,z);
twodfnerr=@(z) pchip(xm,yerr,z);

load G3Cv4up8/G3Cv4up8.mat
n2df=normrnd(pchip(xm,ym,grp.Zfof),pchip(xm,yerr,grp.Zfof));
n2df=max(n2df,pchip(xm,yl(:,1),grp.Zfof));
n2df=min(n2df,pchip(xm,yl(:,2),grp.Zfof));
grp.Mult2dF=grp.VolMult.*n2df;
grp.Mult2dF(grp.Zfof>0.3)=0;
save('G3Cv4up8/G3Cv4up8.mat','grp','-append');
%% galmock environmental properties
for i=1:9
%     mock=load(['201310/mockcat_',num2str(i),'.mat']);
    load(['G3Cv4up8/mockcat_',num2str(i),'.mat']);
%     galmock.IsIterCen=zeros(size(galmock.g));
%     galmock.IsIterCen(grpmock.IterRef)=1;
%     f=(mock.galmock.r_mag<=19.4);
%     galmock.groupID=mock.galmock.groupID(f);
    galmock.EnvLevel=zeros(size(galmock.groupID));
    galmock.EnvLevel(galmock.groupID>0)=1;%0: field; 1: sat; 2: BCG
    %now update BCGRef; done. never do it again!
%     newInd=zeros(size(mock.galmock.g));
%     newInd(f)=1:numel(galmock.g);
%     grpmock.BCGRef=newInd(grpmock.BCGRef);
    if any(grpmock.BCGRef==0)
        error('BCGRef error');
    end
    galmock.EnvLevel(grpmock.BCGRef)=2;
    if numel(galmock.EnvLevel)~=numel(galmock.g)
        error('Env array shape not correct');
    end
    save(['G3Cv4up8/mockcat_',num2str(i),'.mat'],'galmock','grpmock','-append');
end
%% define CentralSampleIter
for i=1:9
    load(['G3Cv4up8/mockcat_',num2str(i),'.mat']);
    galmock.CentralSampleIter=zeros(size(galmock.groupID)); 
    galmock.CentralSampleIter(galmock.IsIterCen>0|galmock.EnvLevel==0)=1;%IterCen or Field
    galmock.CentralSampleIter(galmock.r_mag>19.4)=-1; %undefined.
    save(['G3Cv4up8/mockcat_',num2str(i),'.mat'],'galmock','grpmock','-append');
end
%% define CentralSampleSMpeak
for i=1
    load(['G3Cv4up8/mockcat_',num2str(i),'.mat']);
    galmock.CentralSampleSMpeak=zeros(size(galmock.groupID)); 
    galmock.CentralSampleSMpeak(galmock.EnvLevel==0)=1;%Field
    galmock.CentralSampleSMpeak(grpmock.SMpeakRef)=1;%SMpeak
    galmock.CentralSampleSMpeak(galmock.r_mag>19.4)=-1; %undefined.
    save(['G3Cv4up8/mockcat_',num2str(i),'.mat'],'galmock','grpmock','-append');
end
%% define CentralSampleIterPeak
for i=1
    load(['G3Cv4up8/mockcat_',num2str(i),'.mat']);
    galmock.CentralSampleIterPeak=zeros(size(galmock.groupID)); 
    f=grpmock.IterRef==grpmock.SMpeakRef';
    galmock.CentralSampleIterPeak(grpmock.SMpeakRef(f))=1;%SMpeak
    galmock.CentralSampleIterPeak(galmock.r_mag>19.4)=-1; %undefined.
    save(['G3Cv4up8/mockcat_',num2str(i),'.mat'],'galmock','grpmock','-append');
end
%% define CentralSampleIterNonePeak
for i=1
    load(['G3Cv4up8/mockcat_',num2str(i),'.mat']);
    galmock.CentralSampleIterNonePeak=(galmock.CentralSampleIter>0&galmock.EnvLevel>0&galmock.CentralSampleIterPeak==0); 
    save(['G3Cv4up8/mockcat_',num2str(i),'.mat'],'galmock','grpmock','-append');
end
%% define CentralSampleIter
load(['G3Cv4up8/G3Cv4up8.mat']);
gal.CentralSampleIter=zeros(size(gal.GroupID)); 
gal.CentralSampleIter(gal.IsIterCen>0|gal.EnvLevel==0)=1;%IterCen or Field
gal.CentralSampleIter(gal.Rpetro>19.4)=-1; %undefined.
save(['G3Cv4up8/G3Cv4up8.mat'],'gal','-append');
%% define CentralSampleSMpeak (MMG)
load(['G3Cv4up8/G3Cv4up8.mat']);
gal.CentralSampleSMpeak=zeros(size(gal.GroupID)); 
gal.CentralSampleSMpeak(gal.EnvLevel==0)=1;%Field
gal.CentralSampleSMpeak(grp.SMpeakRef)=1; %group SMpeak
gal.CentralSampleSMpeak(gal.Rpetro>19.4)=-1; %undefined.
save(['G3Cv4up8/G3Cv4up8.mat'],'gal','-append');
%% define CentralSampleIterPeak (MMG and SMpeak)
load(['G3Cv4up8/G3Cv4up8.mat']);
gal.CentralSampleIterPeak=zeros(size(gal.GroupID)); 
f=gal.GalID(grp.SMpeakRef)==grp.IterCenRef;
gal.CentralSampleIterPeak(grp.SMpeakRef(f))=1; %group SMpeak
gal.CentralSampleIterPeak(gal.Rpetro>19.4)=-1; %undefined.
save(['G3Cv4up8/G3Cv4up8.mat'],'gal','-append');
%% define CentralSampleIterNonePeak (MMG and SMpeak)
load(['G3Cv4up8/G3Cv4up8.mat']);
gal.CentralSampleIterNonePeak=(gal.CentralSampleIter>0&gal.EnvLevel>0&gal.CentralSampleIterPeak==0);
save(['G3Cv4up8/G3Cv4up8.mat'],'gal','-append');
%% add SFR
load('G3Cv4up8/G3Cv4up8.mat');
Lines=fits_load_bintable('additional/SFR/EmLinesPhysv04.fits');
%%
Lines.Warning=zeros([numel(Lines.CATAID),1]);
for i=1:numel(Lines.CATAID)
    switch Lines.WARNING_FLAG{i}
        case 'HA    '
            flag=1;
        case 'HB    '
            flag=2;
        case 'Splice'
            flag=3;
        case 'Ok    '
            flag=4;
        case 'XXX   ' %unavailable(can be simply because it's a SDSS galaxy) or problematic
            flag=-1;
        otherwise
            error('unknown warning flag')
    end
    Lines.Warning(i)=flag;
end
%%
Lines.IsAGN=zeros([numel(Lines.CATAID),1]);
for i=1:numel(Lines.CATAID)
    switch Lines.EMLINE_CLASS{i}
        case 'AGN'
            flag=1;
        case 'SF '
            flag=0;
        case 'XXX'
            flag=-1;
        otherwise
            error('unknown agn flag')
    end
    Lines.IsAGN(i)=flag;
end
%%
Lines.FlagBalmer=zeros([numel(Lines.CATAID),1]);
for i=1:numel(Lines.CATAID)
    switch Lines.BALMERDEC_FLAG{i}
        case 'Measured '
            flag=1;
        case 'Estimated'
            flag=0;
        case 'XXX      '
            flag=-1;
        otherwise
            error('unknown balmer flag')
    end
    Lines.FlagBalmer(i)=flag;
end
%%
gal.EmLineInd=zeros([numel(gal.GalID),1]);
for i=1:numel(gal.GalID)
    j=find(gal.GalID(i)==Lines.CATAID);
    if isempty(j)
        j=0;
    end
    if numel(j)>1
        error(['multiple match for id=',num2str(gal.GalID(i))]);
    end
    gal.EmLineInd(i)=j;
end
%%
f=gal.EmLineInd>0;
glist=gal.EmLineInd(f);
gal.EmLineSN=zeros(size(gal.EmLineInd));
gal.FluxHalpha=gal.EmLineSN;
gal.FluxHbeta=gal.EmLineSN;
gal.BalmerDecr=gal.EmLineSN;
gal.FlagBalmer=gal.EmLineSN-1;
gal.FlagWarn=gal.EmLineSN-1;
gal.IsAGN=gal.EmLineSN-1;
gal.LumHalpha=gal.EmLineSN;
gal.SFR=gal.EmLineSN;
%
gal.EmLineSN(f)=Lines.SN_EMI(glist);
gal.FluxHalpha(f)=Lines.F_HALPHA(glist);
gal.FluxHbeta(f)=Lines.F_HBETA(glist);
gal.BalmerDecr(f)=Lines.BALMER_DECREMENT(glist);
gal.FlagBalmer(f)=Lines.FlagBalmer(glist);
gal.LumHalpha(f)=Lines.L_APOBSCOR(glist);
gal.SFR(f)=Lines.SFR(glist);
gal.FlagWarn(f)=Lines.Warning(glist);
gal.IsAGN(f)=Lines.IsAGN(glist);
%
gal.ContSN=zeros(size(gal.EmLineInd));
gal.ContSN(f)=Lines.SN_CONT(glist);
%% change unit:
hubble=0.7;
gal.LumHalpha=gal.LumHalpha*hubble^2; %J/s/h^2
gal.SFR=gal.SFR*1e9*hubble^2; %Msun/h^2/Gyr  (Gyr/h?)
%%
%select non-agn, warning~=HA, SN_LINE>2, -0.01<z<0.31, and has measurement
% Fhalpha>25e-17 erg/cm^2/s, warning~=HB, warning~=Splice.
gal.FlagSFR=gal.EmLineInd>0&gal.FlagWarn~=1&gal.Zspec<0.31&gal.Zspec>-0.01; %good quality SFR data; FlagWarn cut only removes a small amount of gals
gal.FlagSFR=gal.FlagSFR&gal.EmLineSN>2; %high SN, this rejects 40 percent of the above data
gal.FlagSFR=gal.FlagSFR&gal.IsAGN~=1; %not agn
% gal.FlagSFR=gal.FlagSFR&gal.SFR>0; %not necessary. all those negative SFR has already been excluded.
% gal.FlagSFR=gal.FlagSFR&(gal.EnvLevel==0|gal.IsIterCen);
save('G3Cv4up8/G3Cv4up8.mat','gal','-append');
%%
gal.SSFR=gal.SFR./gal.SMsps;
gal=rmfield(gal,{'logSFR','logSSFR'});
save('G3Cv4up8/G3Cv4up8.mat','gal','-append');
%%
for i=1:9
    load(['G3Cv4up8/mockcat_',num2str(i),'.mat']);
    galmock.SSFR=galmock.SFR./galmock.Mstar;
    save(['G3Cv4up8/mockcat_',num2str(i),'.mat'],'galmock','-append');
end
%% add SSFR to groups
ssfr=zeros(max(gal.GalID),1);
tmp=gal.SSFR;
tmp(gal.FlagSFR<1)=nan;
ssfr(gal.GalID)=tmp;
grp.IterCenSSFR=ssfr(grp.IterCenRef);
save('G3Cv4up8/G3Cv4up8.mat','grp','-append');
%% add SSFR to mock groups
for i=1:9
    load(['G3Cv4up8/mockcat_',num2str(i),'.mat']);
    grpmock.SSFRIter=grpmock.SFRIter./grpmock.MstarIter;
    save(['G3Cv4up8/mockcat_',num2str(i),'.mat'],'grpmock','-append');
end
%% add magnitudes to gal
galv6=fits_load_bintable('/work/Projects/Lensing/data/201308/TilingCatv40_kcorr_z00v03.fits');
hubble=1.0; %assume H=100, so all the magnitudes and luminosities should retain h in the units.
dm=@(z) 5*log10(lum_dist(0.3,z)*1e6/hubble/10); %should use asinh() DM??
kcorr=@(p1,p2,p3,p4,p5,z) p1.*z.^4+p2.*z.^3+p3.*z.^2+p4.*z+p5;
zgal=galv6.Z_TONRY;
DM=dm(zgal);DM(zgal<=0)=0;
Umodel_abs=galv6.U_MODEL-DM-kcorr(galv6.PCOEFF_U_1,galv6.PCOEFF_U_2,galv6.PCOEFF_U_3,galv6.PCOEFF_U_4,galv6.PCOEFF_U_5,zgal);
Gmodel_abs=galv6.G_MODEL-DM-kcorr(galv6.PCOEFF_G_1,galv6.PCOEFF_G_2,galv6.PCOEFF_G_3,galv6.PCOEFF_G_4,galv6.PCOEFF_G_5,zgal);
Rmodel_abs=galv6.R_MODEL-DM-kcorr(galv6.PCOEFF_R_1,galv6.PCOEFF_R_2,galv6.PCOEFF_R_3,galv6.PCOEFF_R_4,galv6.PCOEFF_R_5,zgal);
Imodel_abs=galv6.I_MODEL-DM-kcorr(galv6.PCOEFF_I_1,galv6.PCOEFF_I_2,galv6.PCOEFF_I_3,galv6.PCOEFF_I_4,galv6.PCOEFF_I_5,zgal);
Zmodel_abs=galv6.Z_MODEL-DM-kcorr(galv6.PCOEFF_Z_1,galv6.PCOEFF_Z_2,galv6.PCOEFF_Z_3,galv6.PCOEFF_Z_4,galv6.PCOEFF_Z_5,zgal);
gal.TilingInd=NaN(size(gal.GalID));
for i=1:numel(gal.GalID)
    ind=find(galv6.CATAID==gal.GalID(i));
    if numel(ind)==1
        gal.TilingInd(i)=ind;
    end
end
flag=~isnan(gal.TilingInd);
arr_nan=NaN(size(gal.GalID));
gal.Umodel_abs=arr_nan;
gal.Umodel_abs(flag)=Umodel_abs(gal.TilingInd(flag));
gal.Gmodel_abs=arr_nan;
gal.Gmodel_abs(flag)=Gmodel_abs(gal.TilingInd(flag));
gal.Rmodel_abs=arr_nan;
gal.Rmodel_abs(flag)=Rmodel_abs(gal.TilingInd(flag));
gal.Imodel_abs=arr_nan;
gal.Imodel_abs(flag)=Imodel_abs(gal.TilingInd(flag));
gal.Zmodel_abs=arr_nan;
gal.Zmodel_abs(flag)=Zmodel_abs(gal.TilingInd(flag));
save('G3Cv4up8/G3Cv4up8.mat','gal','-append');
%%
clear
load G3Cv4up8/G3Cv4up8.mat
%gal.Umodel_abs(gal.Umodel_abs<-100)=NaN;
f=gal.Umodel_abs>-10|gal.Gmodel_abs>-10|gal.Rmodel_abs>-10|gal.Imodel_abs>-10|gal.Zmodel_abs>-10; %this remove ~100 gals.
gal.Umodel_abs(f)=nan;
gal.Gmodel_abs(f)=nan;
gal.Rmodel_abs(f)=nan;
gal.Imodel_abs(f)=nan;
gal.Zmodel_abs(f)=nan;
save('G3Cv4up8/G3Cv4up8.mat','gal','-append');
%% add magnitudes to mock
for i=1:9
    mock=load(['201310/mockcat_',num2str(i),'.mat']);
    load(['G3Cv4up8/mockcat_',num2str(i),'.mat']);
    galmock.OldGalInd=find(mock.galmock.r_mag<=19.4);
    if numel(galmock.OldGalInd)~=numel(galmock.galID)
        error('OldGalInd err');
    end
    galmock.Umodel_abs=mock.galmock.Umag(galmock.OldGalInd);
    galmock.Gmodel_abs=mock.galmock.Gmag(galmock.OldGalInd);
    galmock.Rmodel_abs=mock.galmock.Rmag(galmock.OldGalInd);
    galmock.Imodel_abs=mock.galmock.Imag(galmock.OldGalInd);
    galmock.Zmodel_abs=mock.galmock.Zmag(galmock.OldGalInd);
    save(['G3Cv4up8/mockcat_',num2str(i),'.mat'],'galmock','-append');
end
%% save for Wenting
clear;
clc;
load G3Cv4up8/G3Cv4up8.mat
f=gal.CentralSampleIter>0&~isnan(gal.SMsps);
x=[gal.RA(f),gal.DEC(f),gal.Rpetro(f), gal.Zspec(f),gal.SMsps(f),gal.FlagSFR(f)>0&gal.SSFR(f)>10^-1.5];
save GAMAcentral.dat x -ascii
%% random sample (random positions only)
clear;
clc;
load G3Cv4up8/G3Cv4up8.mat
f=gal.CentralSampleIter>0&~isnan(gal.SMsps);
data1.RA=gal.RA(f);
data1.DEC=gal.DEC(f);
[g1,g2,g3]=split_mockcat(data1);
rnd=GAMArand([numel(g1.RA),numel(g2.RA),numel(g3.RA)],-10);
rand1.RA=[rnd{1}(:,1);rnd{2}(:,1);rnd{3}(:,1)];
rand1.DEC=[rnd{1}(:,2);rnd{2}(:,2);rnd{3}(:,2)];
x=[rand1.RA,rand1.DEC,gal.Rpetro(f), gal.Zspec(f),gal.SMsps(f),gal.FlagSFR(f)>0&gal.SSFR(f)>10^-1.5];
save GAMAcentral_random.dat x -ascii
%% Update sample flags
%FlagSFR: select non-agn, warning~=HA, SN_LINE>2, -0.01<z<0.31, and has match
%hierarchical FlagSFR: >=1: generally good; >=2: cleaner sample; >=3: super
load G3Cv4up8/G3Cv4up8.mat
f=gal.FlagSFR>0&gal.FlagWarn==4;%restrain to FlagWarn=OK.
f=f&gal.EmLineSN>5&~isnan(gal.SMsps); %higher SN cut, and have SM value
gal.FlagSFR=int32(gal.FlagSFR);
gal.FlagSFR(f)=2;
f=f&gal.FlagBalmer==1; %only measured Balmer line (not estimated; i.e, with both Ha and Hb measurement).
f=f&gal.IsAGN==0; %SFgalaxy. almost all the Balmer measured galaxies are SFgalaxies. only a tiny fraction is AGN.
gal.FlagSFR(f)=3; 
gal.ActiveSample=zeros(size(gal.FlagSFR))-1; %default: -1; not inside the samples 
gal.ActiveSample(gal.CentralSampleIter>0&gal.FlagSFR>=2&gal.SSFR>10^-1.5)=1; %active ones
gal.ActiveSample(gal.CentralSampleIter>0&gal.FlagSFR>=2&gal.SSFR<10^-1.5)=0; %passive ones
save('G3Cv4up8/G3Cv4up8.mat','gal','-append');
%% FluxHalpha>25
%select non-agn, warning==OK or XXX, -0.01<z<0.31, and has measurement in
%SFR and SM
% Fhalpha>25e-17 erg/cm^2/s
f=gal.EmLineInd>0&(gal.FlagWarn==4|gal.FlagWarn==-1)&gal.Zspec<0.31&gal.Zspec>-0.01; %good quality SFR data; FlagWarn cut removes 10% of galaxies;
f=f&gal.FluxHalpha>25; %flux cut, removes 35% from above
f=f&gal.IsAGN~=1; %not agn, removes 5% from above
f=f&~isnan(gal.SMsps); %have SM value, removes 4% from above
gal.FlagSFRFluxCut=f;
gal.ActiveFullSampleFluxCut=zeros(size(f))-1; %default: -1; not inside the samples 
gal.ActiveFullSampleFluxCut(f>0&gal.SSFR>10^-1.5)=1; %active ones
gal.ActiveFullSampleFluxCut(f>0&gal.SSFR<10^-1.5)=0; %passive ones
save('G3Cv4up8/G3Cv4up8.mat','gal','-append');
%% SN_CONT>3 cut
%select non-agn, warning==OK or XXX, -0.01<z<0.31, and has measurement in
%SFR and SM
% ContSN>3
f=gal.EmLineInd>0&(gal.FlagWarn==4|gal.FlagWarn==-1)&gal.Zspec<0.31&gal.Zspec>-0.01; %good quality SFR data; FlagWarn cut removes 10% of galaxies;
f=f&gal.ContSN>3; %SN_Continuum cut, removes 25% from above
f=f&gal.IsAGN~=1; %not agn, removes 8% from above
f=f&~isnan(gal.SMsps); %have SM value, removes 3% from above
gal.FlagSFRContCut=f;
gal.ActiveFullSampleSNContCut=zeros(size(f))-1; %default: -1; not inside the samples 
gal.ActiveFullSampleSNContCut(f>0&gal.SSFR>10^-1.5)=1; %active ones
gal.ActiveFullSampleSNContCut(f>0&gal.SSFR<10^-1.5)=0; %passive ones
save('G3Cv4up8/G3Cv4up8.mat','gal','-append');
%% Combined SN: SN_EMI>3|(SN_Cont>3&SN_EMI<0), almost the same central sample as ContCut
%select non-agn, warning==OK or XXX, -0.01<z<0.31, and has measurement in
%SFR and SM
% ContSN>3
f=gal.EmLineInd>0&(gal.FlagWarn==4|gal.FlagWarn==-1)&gal.Zspec<0.31&gal.Zspec>-0.01; %good quality SFR data; FlagWarn cut removes 10% of galaxies;
f=f&(gal.EmLineSN>3|(gal.EmLineSN<0&gal.ContSN>3)); %SN_EMI>3 or Cont>3 when no EMI_SN availabel, removes 30% from above
f=f&gal.IsAGN~=1; %not agn, removes 8% from above
f=f&~isnan(gal.SMsps); %have SM value, removes 3% from above
gal.FlagSFRCombCut=f;
gal.ActiveFullSampleSNCombCut=zeros(size(f))-1; %default: -1; not inside the samples 
gal.ActiveFullSampleSNCombCut(f>0&gal.SSFR>10^-1.5)=1; %active ones
gal.ActiveFullSampleSNCombCut(f>0&gal.SSFR<10^-1.5)=0; %passive ones
save('G3Cv4up8/G3Cv4up8.mat','gal','-append');
%% Halpha Lum Cut: LumHalpha>1e34
%select non-agn, warning==OK or XXX, -0.01<z<0.31, and has measurement in
%SFR and SM
% ContSN>3
f=gal.EmLineInd>0&(gal.FlagWarn==4|gal.FlagWarn==-1)&gal.Zspec<0.31&gal.Zspec>-0.01; %good quality SFR data; FlagWarn cut removes 10% of galaxies;
f=f&gal.IsAGN~=1; %not agn, removes 8% from above
f=f&~isnan(gal.SMsps); %have SM value, removes 3% from above
% f=f&gal.SFR>0; %have SFR value, removes 30% from above
gal.FlagSFRBaseCut=f;
gal.HaLumSample=zeros(size(f))-1; %default: -1; not inside the samples 
gal.HaLumSample(f>0&gal.LumHalpha>1e34)=1; %active ones
gal.HaLumSample(f>0&gal.LumHalpha<1e34)=0; %passive ones
save('G3Cv4up8/G3Cv4up8.mat','gal','-append');
%% add alternative VelDisp
load G3Cv4up8/G3Cv4up8.mat
grp.VelDisp2=sqrt(grp.VelDispRaw.^2.*(grp.Mult-1)./grp.Mult-grp.VelErr.^2);
grp.VelDisp3=sqrt(grp.VelDispRaw.^2.*grp.Mult./(grp.Mult-1)-grp.VelErr.^2);
grp.VelDisp2=real(grp.VelDisp2);
grp.VelDisp3=real(grp.VelDisp3);
save('G3Cv4up8/G3Cv4up8.mat','grp','-append');
for i=1:9
    load(['G3Cv4up8/mockcat_',num2str(i),'.mat']);
    grpmock.VelDisp2=sqrt(grpmock.VelDispRaw.^2.*(grpmock.Mult-1)./grpmock.Mult-grpmock.VelErr.^2);
    grpmock.VelDisp3=sqrt(grpmock.VelDispRaw.^2.*grpmock.Mult./(grpmock.Mult-1)-grpmock.VelErr.^2);
    grpmock.VelDisp2=real(grpmock.VelDisp2);
    grpmock.VelDisp3=real(grpmock.VelDisp3);
    save(['G3Cv4up8/mockcat_',num2str(i),'.mat'],'grpmock','-append');
end