clc;
clear;
clear global;

disp('loading GAMA...');
cd /work/Projects/Lensing/data
load G3Cv4up8/G3Cv4up8.mat
%%
grp.GroupID=int32(grp.GroupID);
grp.Mult=int32(grp.Mult);
grp.Mult2dF=int32(grp.Mult2dF);
grp.IterCenRef=int32(grp.IterCenRef);
grp.BCGRef=int32(grp.BCGRef);
grp.LinkTot=int32(grp.LinkTot);
[grp09,grp12,grp15]=split_mockcat(grp);
gr={grp09,grp12,grp15};
%%
names=fieldnames(gr{1});
ncol=numel(names);
hdf5write('groupcatv4up8.hdf5','/info','gama-I G3Cv4 updated for lensingv8.0 groups for 3 sky patches');
for sky=1:3
for i=1:ncol
hdf5write('groupcatv4up8.hdf5',['/sky',num2str(sky-1),'/',names{i}],gr{sky}.(names{i}),'writemode','append')
end
end
%% gal the same as G3Cv4up6
% !ln GAMAgalv4up6.hdf5 GAMAgalv4up8.hdf5 -s
gal.EnvLevel=int32(gal.EnvLevel); %environment: 0: field; 1: sat; 2: bcg.
gal.IsIterCen=int32(gal.IsIterCen); % 0: not itercen; >0: GroupID in which it is itercen
gal.GAMAID=int32(gal.GalID);
gal.Source=int32(gal.S);
gal.GroupID=int32(gal.GroupID);
gal.FlagBalmer=int32(gal.FlagBalmer);
gal.FlagWarn=int32(gal.FlagWarn);
gal.FlagSFR=int32(gal.FlagSFR);
gal.IsAGN=int32(gal.IsAGN);
gal.EmLineInd=int32(gal.EmLineInd);
gal.CentralSampleIter=int32(gal.CentralSampleIter);
gal.ActiveSample=int32(gal.ActiveSample);
gal.FlagSFRFluxCut=int32(gal.FlagSFRFluxCut);
gal.FlagSFRContCut=int32(gal.FlagSFRContCut);
gal.FlagSFRCombCut=int32(gal.FlagSFRCombCut);
gal.FlagSFRBaseCut=int32(gal.FlagSFRBaseCut);
gal.ActiveFullSampleFluxCut=int32(gal.ActiveFullSampleFluxCut);
gal.ActiveFullSampleSNContCut=int32(gal.ActiveFullSampleSNContCut);
gal.ActiveFullSampleSNCombCut=int32(gal.ActiveFullSampleSNCombCut);
gal.HaLumSample=int32(gal.HaLumSample);
% gal.BinID=int32(gal.BinID);
gal.CentralSampleSMpeak=int32(gal.CentralSampleSMpeak);
gal.CentralSampleIterPeak=int32(gal.CentralSampleIterPeak);
gal.CentralSampleIterNonePeak=int32(gal.CentralSampleIterNonePeak);
gal=rmfield(gal,{'S','GalID'});
[g{1},g{2},g{3}]=split_mockcat(gal);
%%
names=fieldnames(g{1});
ncol=numel(names);
hdf5write('GAMAgalv4up8.hdf5','/info','gama galaxies G3Cv4 for 3 sky patches, updated StellarMass from G3Cv6, updated SFR from EmLinePhys');
for sky=1:3
for i=1:ncol
hdf5write('GAMAgalv4up8.hdf5',['/sky',num2str(sky-1),'/',names{i}],g{sky}.(names{i}),'writemode','append')
end
end
%%
disp('loading shear...');
% e{1}=loadShearTBL(1,[1,2,40,41,42,53,57,58]);% ra,dec,chi1,chi2,sigma_e,z,z-sig,z+sig
% e{2}=loadShearTBL(2,[1,2,40,41,42,53,57,58]);
% e{3}=loadShearTBL(3,[1,2,40,41,42,53,57,58]);
% save shearmap.mat e

load shearmap
for sky=1:3
    e{sky}(:,5)=e{sky}(:,5)*2;
end  %change the error to error in chi (or e) rather than in epsilon


e1=e{1};
e2=e{2};
e3=e{3};
save shearmap_ZEBRA e1 e2 e3
%%
load GAMA_SRC_matchZ
for i=1:3
    e{i}(source_idmatch{i},6:8)=[gama_zmatch{i}(:,1),gama_zmatch{i}(:,1)-gama_zmatch{i}(:,2)/3e5,gama_zmatch{i}(:,1)+gama_zmatch{i}(:,2)/3e5];
end
e1=e{1};
e2=e{2};
e3=e{3};
save shearmap_ZEBRAGAMA e1 e2 e3
%% to do the test for Peacock
xbrd=[1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11];
binID=zeros(size(grp.GroupID)); %default in bin 0
for i=1:size(xbrd,1)
    f=grp.Mult>2&grp.IterCenSMsps>xbrd(i,1)&grp.IterCenSMsps<=xbrd(i,2);
    binID(f)=i;
end
gal.BinID=zeros(size(gal.GalID))-1; %undefined, not a central.
gal.BinID(gal.CentralSampleIter>0)=0; %bin 0, inside CentralSampleIter, but not in other bins
for i=1:numel(grp.GroupID)
    gal.BinID(gal.GalID==grp.IterCenRef(i))=binID(i);
end
