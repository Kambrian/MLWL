for mockid=1:9
%%
cd /home/kam/Projects/Lensing/data
grpmock=loadGAMAcsv(['G3Cv4/mock/G3CMock',num2str(mockid),'FoFGroup194v04.dat'],4);
halomock=loadGAMAcsv(['G3Cv4/mock/G3CMock',num2str(mockid),'HaloGroup194v04.dat'],4);
galmock=loadGAMAcsv(['G3Cv4/mock/G3CMock',num2str(mockid),'Galv04.dat'],4);
refmock=loadGAMAcsv(['G3Cv4/mock/G3CMock',num2str(mockid),'Ref194v04.dat'],4);
Ac=-1.2;An=20.7;Az=2.3;
Bc=0.94;Bn=-0.67;Bz=0.16;
% halomock.MassProxy=halomock.MassProxy.*(Ac+An./sqrt(halomock.Mult)+Az./sqrt(halomock.Zfof));
% halomock.MassProxy=halomock.MassProxy.*dynmass_calibration(halomock);   %hardly any difference
% halomock.TotFluxProxy=halomock.TotFluxProxy.*(Bc+Bn./sqrt(halomock.Mult)+Bz./sqrt(halomock.Zfof));
% [halomock.LumMass,C1,A1]=luminosity_mass(halomock);
%[grpmock.DynMass,CD,AD]=dynamical_mass(grpmock,mockmass.MIter);
% grpmock.MassProxy=grpmock.MassProxy.*(Ac+An./sqrt(grpmock.Mult)+Az./sqrt(grpmock.Zfof));
% grpmock.MassProxy=grpmock.MassProxy.*dynmass_calibration(grpmock);
% grpmock.TotFluxProxy=grpmock.TotFluxProxy.*(Bc+Bn./sqrt(grpmock.Mult)+Bz./sqrt(grpmock.Zfof));
% [grpmock.LumMass,C2,A2]=luminosity_mass(grpmock);

% grp=loadGAMAcsv('G3Cv4/group/G3CFoFGroup194v04.dat',4);
% ref=loadGAMAcsv('G3Cv4/group/G3CRef194v04.dat',4);
% grp.MassProxy=grp.MassProxy.*(Ac+An./sqrt(grp.Mult)+Az./sqrt(grp.Zfof));
% grp.TotFluxProxy=grp.TotFluxProxy.*(Bc+Bn./sqrt(grp.Mult)+Bz./sqrt(grp.Zfof));
% [grp.LumMass,C3,A3]=luminosity_mass(grp);
% grp.HybMass=sqrt(grp.MassProxy.*grp.LumMass);
% 
% gal=loadGAMAcsv('G3Cv4/group/G3CGalv04.dat',4);

linkmock=[refmock.GalID,refmock.GroupID];  %galID: index (row number) for galaxies in the galaxy file
linkmock=sortrows(linkmock,2);
% link=[ref.GalID,ref.GroupID];
% link=sortrows(link,2);
%%
ngal=numel(galmock.GalID);
ngrp=numel(grpmock.GroupID);
%% grp2halo
GrpOffset=sum(~linkmock(:,2));
grp2halo=zeros(ngrp,1);
for grpid=1:ngrp
    gallist=linkmock(GrpOffset+(1:grpmock.Mult(grpid)),1);
    halolist=galmock.HaloID(gallist);
    [haloid,nshare]=mode(halolist);
    if nshare>grpmock.Mult(grpid)/2
        grp2halo(grpid)=haloid;
    end
    GrpOffset=GrpOffset+grpmock.Mult(grpid);
end
%% halo2grp
linkmockhalo=[galmock.GalID,galmock.HaloID];
linkmockhalo=sortrows(linkmockhalo,2);
nhalo=max(galmock.HaloID);
halo2grp=zeros(nhalo,1);
haloid=linkmockhalo(1,2);
gallist=linkmockhalo(1,1);
for i=2:ngal
    if linkmockhalo(i,2)==haloid
        gallist=[gallist;linkmockhalo(i,1)];
        continue;
    end
    %finish previous halo
    grplist=refmock.GroupID(gallist);
    [grpid,nshare]=mode(grplist);
    if nshare>numel(gallist)/2
        halo2grp(haloid)=grpid;
    end
    %start a new halo
    gallist=linkmockhalo(i,1);
    haloid=linkmockhalo(i,2);
end
%finish last halo
grplist=refmock.GroupID(gallist);
[grpid,nshare]=mode(grplist);
if nshare>numel(gallist)
    halo2grp(haloid)=grpid;
end
%%
grpmass=zeros(ngrp,1);
for grpid=1:ngrp
    haloid=grp2halo(grpid);
    if haloid>0 %found a halo
        if halo2grp(haloid)==grpid
            grpmass(grpid)=halomock.HaloMass(halomock.HaloID==haloid);
        end
    end
end
%%
% cd /home/kam/Projects/Lensing/data
% grpmock=loadGAMAcsv('G3Cv4/mock/G3CMock1FoFGroup194v04.dat',4);
[grpmock.DynMass,CD1{mockid},CD2{mockid},AD1{i},AD2{i}]=dynamical_mass(grpmock,grpmass);
A1{mockid}=median(grpmass(grpmass>0&grpmock.Mult>4)./grpmock.MassProxy(grpmass>0&grpmock.Mult>4));
A2{mockid}=median(grpmass(grpmass>0&grpmock.Mult>4))./median(grpmock.MassProxy(grpmass>0&grpmock.Mult>4));
end
%%
CD1=cell2mat(reshape(CD1,1,1,9));
CD2=cell2mat(reshape(CD2,1,1,9));
C=[19.0,10.8,12.0,12.6;
19.5,10.5,11.1,10.4;
21.5,10.3,8.6,8.3;
17.4,6.1,5.4,5.6;
];
mean(CD1,3)
mean(CD2,3)