clc;
clear;
clear global;

disp('loading GAMA...');
cd /work/Projects/Lensing/data
load G3Cv6.mat

grp.GroupID=int32(grp.GroupID);
grp.Mult=int32(grp.Mult);
grp.IterCenRef=int32(grp.IterCenRef);
grp.BCGRef=int32(grp.BCGRef);
grp.LinkTot=int32(grp.LinkTot);
[grp09,grp12,grp15]=split_mockcat(grp);
gr={grp09,grp12,grp15};
%%
names=fieldnames(gr{1});
ncol=numel(names);
hdf5write('groupcatv6.hdf5','/info','gama-II G3Cv6 groups for 3 sky patches');
for sky=1:3
for i=1:ncol
hdf5write('groupcatv6.hdf5',['/sky',num2str(sky-1),'/',names{i}],gr{sky}.(names{i}),'writemode','append')
end
end
%%
gal.EnvLevel=int32(gal.EnvLevel); %environment: 0: field; 1: sat; 2: bcg.
gal.IsIterCen=int32(gal.IsIterCen); % 0: not itercen; >0: GroupID in which it is itercen
gal.IsBCG=int32(gal.IsBCG); 
gal.GAMAID=int32(gal.CATAID);
gal.Source=int32(gal.SURVEY_CODE);
gal.GroupID=int32(gal.GroupID);
gal=rmfield(gal,{'SURVEY_CODE','CATAID'});
[g{1},g{2},g{3}]=split_mockcat(gal);
%%
names=fieldnames(g{1});
ncol=numel(names);
hdf5write('GAMAgalv6.hdf5','/info','gama-II galaxies G3Cv6 for 3 sky patches');
for sky=1:3
for i=1:ncol
hdf5write('GAMAgalv6.hdf5',['/sky',num2str(sky-1),'/',names{i}],g{sky}.(names{i}),'writemode','append')
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