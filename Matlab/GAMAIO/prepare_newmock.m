%% to collect and put the information together for the new mocks
MagSun_r=4.67;
cd /work/Projects/Lensing/data
Ac=2.0;An=17.9;Az=1.5; %r<19.8
% Bc=0.65;Bn=-0.50;Bz=0.22;%r<19.8
Bc=0.086;Bn=-0.37;Bz=0.35; %new calibration, r<19.8
%%
for i=1:9
    disp(['processing mock ',num2str(i),'...']);
grpmockall{i}=loadGAMAcsv(['201310/newmocks/FoFtab',num2str(i),'.csv'],6);
% halomockall{i}=loadGAMAcsv('201310/mock/G3CMock1HaloGroup194v04.dat',4);
galmockall{i}=loadGAMAcsv(['201310/galmock_',num2str(i-1),'.dat'],6);
galmockall{i}.Mcolor=10.^(1.15+0.70*(galmockall{i}.g-galmockall{i}.i)-0.4*galmockall{i}.i)*0.7^2;%Msun/h^2, what Taylor estimates
refmockall{i}=loadGAMAcsv(['201310/newmocks/FoFrefs',num2str(i),'.csv'],6);
grpmockall{i}.TotFluxProxyRaw=10.^(-0.4*(grpmockall{i}.TotFluxInt-MagSun_r));%Lsun/h^2
grpmockall{i}=rmfield(grpmockall{i},'TotFluxInt');
grpmockall{i}=mvfield(grpmockall{i},{'Mass','MedianZ'},{'MassProxyRaw','Zfof'});
grpmockall{i}.DynMass=grpmockall{i}.MassProxyRaw.*(Ac+An./sqrt(grpmockall{i}.Mult)+Az./sqrt(grpmockall{i}.Zfof));
grpmockall{i}.Luminosity=grpmockall{i}.TotFluxProxyRaw.*(Bc+Bn./sqrt(grpmockall{i}.Mult)+Bz./sqrt(grpmockall{i}.Zfof));
[grpmockall{i}.LumMass,C1,A1]=luminosity_mass(grpmockall{i});

linkmockall{i}=[refmockall{i}.galID,refmockall{i}.groupID];  %galID: index (row number) for galaxies in the galaxy file
linkmockall{i}=sortrows(linkmockall{i},2);

for j=1:numel(grpmockall{i}.Gnum)
    %same galID, and within 10arcsec of ra,dec
    ind=find(grpmockall{i}.IterRef(j)==galmockall{i}.galID&abs(grpmockall{i}.IterCenRA(j)-galmockall{i}.RA)<5e-5&abs(grpmockall{i}.IterCenDEC(j)-galmockall{i}.DEC)<5e-5);
    if numel(ind)~=1
        error(['Failed to match Iter central gal, for grp',num2str(j),', matched gals:',num2str(ind)]);
    end
    grpmockall{i}.IterRef(j)=ind;
    
    ind=find(grpmockall{i}.BCGRef(j)==galmockall{i}.galID&abs(grpmockall{i}.BCGRA(j)-galmockall{i}.RA)<5e-5&abs(grpmockall{i}.BCGDEC(j)-galmockall{i}.DEC)<5e-5);
    if numel(ind)~=1
        error(['Failed to match BCG central gal, for grp',num2str(j),', matched gals:',num2str(ind)]);
    end
    grpmockall{i}.BCGRef(j)=ind;
end
grpmockall{i}.MIter=galmockall{i}.Mhalo(grpmockall{i}.IterRef);
grpmockall{i}.McolorIter=galmockall{i}.Mcolor(grpmockall{i}.IterRef);
grpmockall{i}.MstarIter=galmockall{i}.Mstar(grpmockall{i}.IterRef);
grpmock=grpmockall{i};
refmock=refmockall{i};
galmock=galmockall{i};
linkmock=linkmockall{i};
save(['201310/mockcat_',num2str(i),'.mat'],'grpmock','refmock','galmock','linkmock');
end
%% to fix the bug in TotFluxInt to TotFluxRaw 
% for i=1:9
%     load(['201310/mockcat_',num2str(i),'.mat']);
%     grpmock.TotFluxProxyRaw=10.^(-0.4*(grpmock.TotFluxProxyRaw-MagSun_r));%Lsun/h^2,  
%     grpmock.Luminosity=grpmock.TotFluxProxyRaw.*(Bc+Bn./sqrt(grpmock.Mult)+Bz./sqrt(grpmock.Zfof));
%     [grpmock.LumMass,C1,A1]=luminosity_mass(grpmock);
%     save(['201310/mockcat_',num2str(i),'.mat'],'grpmock','-append');
% end
%% change z to Z
for i=1:9
    load(['201310/mockcat_',num2str(i),'.mat']);
    galmock=mvfield(galmock,'z','Z');
    save(['201310/mockcat_',num2str(i),'.mat'],'galmock','-append');
end
%% (save for Lingyu)
for i=1:9
    load(['201310/mockcat_',num2str(i),'.mat']);
    fp=fopen(['201310/Lingyu/grpmore_',num2str(i),'.dat'],'w');
    fprintf(fp,'Gnum,IterIndex,BCGIndex,DynMass_Msunh,Luminosity_Lsunh2,LumMass_Msunh\n');
    for j=1:numel(grpmock.Gnum)
        fprintf(fp,'%d,%d,%d,%g,%g,%g\n',grpmock.Gnum(j),grpmock.IterRef(j),grpmock.BCGRef(j),grpmock.DynMass(j),grpmock.Luminosity(j),grpmock.LumMass(j));
    end
    fclose(fp);
end
%% Append VolMult, and new lum calibration
% Nav=load('GAMAII-Ndens-z.txt','-ascii');
% avNdens=@(z) pchip(Nav(:,1),Nav(:,2),z);

% declim=[-11.25,-2.75;-2,4;-3,3;-3,3];dra=[11.5;14;14;14];Area=pi/180*sum((sind(declim(:,2))-sind(declim(:,1))).*dra);%much larger area!!
% zbin=[0:0.03:0.52]';
% zmid=(zbin(1:end-1)+zbin(2:end))/2;
% rbin=comoving_dist(0.3,0.7,zbin);
% Vbin=Area/3*diff(rbin.^3);

Bc=0.086;Bn=-0.37;Bz=0.35; %new calibration, r<19.8
for i=1:9
    load(['201310/mockcat_',num2str(i),'.mat']);
%     Nmock=histc(galmock.Z,zbin);
%     nmock=Nmock(1:end-1)./Vbin; %comoving density
%     grpmock.VolMult=grpmock.Mult./pchip(zmid,nmock,grpmock.Zfof);% calibrate with its own dN/dV
    
%     grpmock.VolMult=grpmock.Mult./avNdens(grpmock.Zfof);
    
    grpmock.Luminosity=grpmock.TotFluxProxyRaw.*(Bc+Bn./sqrt(grpmock.Mult)+Bz./sqrt(grpmock.Zfof)); %redo the calibration
    [grpmock.LumMass,C1,A1]=luminosity_mass(grpmock);
    save(['201310/mockcat_',num2str(i),'.mat'],'grpmock','-append');
end
%% load the new mocks
for i=1:9
    mock{i}=load(['201310/mockcat_',num2str(i),'.mat']);
end
%% add SFR
for i=1:9
    disp(['processing mock ',num2str(i),'...']);
    load(['201310/mockcat_',num2str(i),'.mat']);
    galmocktmp=loadGAMAcsv(['201310/galmock_',num2str(i-1),'.dat'],6);
    galmock.SFR=galmocktmp.SFR;
    grpmock.SFRIter=galmock.SFR(grpmock.IterRef);
    save(['201310/mockcat_',num2str(i),'.mat'],'grpmock','galmock','-append');
end
%% redo refmock and linkmock: 26-Nov-2013
for i=1:9
    disp(['processing mock ',num2str(i),'...']);
    load(['201310/mockcat_',num2str(i),'.mat']);
    refmock=loadGAMAcsv(['201310/FoFRef/NewRef_',num2str(i-1),'.dat'],6);
    linkmock=[refmock.galInd,refmock.groupID];  %galInd: index (row number) for galaxies in the galaxy file, starting from 1
    linkmock=sortrows(linkmock,2);
    hostmock=loadGAMAcsv(['201310/FoFRef/GalHosts_',num2str(i-1),'.dat'],6);
    galmock.groupID=hostmock.groupID;
    save(['201310/mockcat_',num2str(i),'.mat'],'galmock','refmock','linkmock','-append');
end
%% add more magnitudes
for i=1:9
    disp(['processing mock ',num2str(i),'...']);
    load(['201310/mockcat_',num2str(i),'.mat']);
    galmocktmp=loadGAMAcsv(['201310/galmock_',num2str(i-1),'.dat'],6);
    galmock.Umag=galmocktmp.Umag;
    galmock.Gmag=galmocktmp.Gmag;
    galmock.Rmag=galmocktmp.Rmag;
    galmock.Imag=galmocktmp.Imag;
    galmock.Zmag=galmocktmp.Zmag;
    save(['201310/mockcat_',num2str(i),'.mat'],'galmock','-append');
end