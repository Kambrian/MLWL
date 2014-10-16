% GAMA WL initialization script
%
%
% clc
clear global
global macros

if exist('macro','var')
    macros=macro;
else
    macros=default_WLparam();
end
%%
cd /work/Projects/Lensing/data

global grp e
disp('loading GAMA...');

switch macros.LENSCATID
    case 1
grp=cell(1,3);
grpmock=loadGAMAcsv('G3Cv4/mock/G3CMock1FoFGroup194v04.dat',4);
[grp{1},grp{2},grp{3}]=split_mockcat(grpmock);
    otherwise 
 Ac=-1.2;An=20.7;Az=2.3;
Bc=0.94;Bn=-0.67;Bz=0.16;       
grp=loadGAMAcsv('G3Cv4/group/G3CFoFGroup194v04.dat',4);
grp.MassProxy=grp.MassProxy.*(Ac+An./sqrt(grp.Mult)+Az./sqrt(grp.Zfof));
grp.TotFluxProxy=grp.TotFluxProxy.*(Bc+Bn./sqrt(grp.Mult)+Bz./sqrt(grp.Zfof));
[grp.LumMass,C3,A3]=luminosity_mass(grp);
grp.HybMass=sqrt(grp.MassProxy.*grp.LumMass);
%cattmp=loadGAMAcsv('G3Cv4/group/G3CFoFGroup194v04.dat',4);
cattmp=grp;
grp=[];
[grp{1},grp{2},grp{3}]=split_mockcat(cattmp);
% g{1}=loadGAMAcsv('groupsv2/cutcatG09.csv');

% g{2}=loadGAMAcsv('groupsv2/cutcatG12.csv');

% g{3}=loadGAMAcsv('groupsv2/cutcatG15.csv');
% cat=grp;
% cat=structcat([cat09;cat12;cat15]);
% gal=structcat([g09;g12;g15]);
end

if isfield(macros,'RTRIM')
grp=gama_trim(grp,macros.RTRIM);
end

if macros.LENSCATID<0 %random
    grpcenters=GAMArand([numel(grp{1}.Gnum),numel(grp{2}.Gnum),numel(grp{3}.Gnum)],macros.LENSCATID);
else %real
    grpcenters=cell(3,1);
    switch macros.BCG_CENTER
        case 0
            for sky=1:3
                grpcenters{sky}=[grp{sky}.CenRA,grp{sky}.CenDEC];
            end
        case 1
            for sky=1:3
                grpcenters{sky}=[grp{sky}.BCGRA,grp{sky}.BCGDEC];
            end
        case 2
            for sky=1:3
                grpcenters{sky}=[grp{sky}.IterCenRA,grp{sky}.IterCenDEC];
            end
        otherwise
            error('unknown center type');
    end
    if macros.LENSCATID>1 %perm
        permarr=GAMArand([numel(grp{1}.Gnum),numel(grp{2}.Gnum),numel(grp{3}.Gnum)],macros.LENSCATID);
        for sky=1:3
            grpcenters{sky}=grpcenters{sky}(permarr);
        end
    end
end
% if isfield(macros,'BCG_Z')
% if macros.BCG_Z
%     g=cell(1,3);
%     g{1}=loadGAMAcsv('groupsv2/cutcatG09.csv');
%     g{2}=loadGAMAcsv('groupsv2/cutcatG12.csv');
%     g{3}=loadGAMAcsv('groupsv2/cutcatG15.csv');
%     for sky=1:3
%         cat{sky}.BCGz=cat{sky}.MedianZ;
%         for gid=1:numel(cat{sky}.Gnum)
%             cat{sky}.BCGz(gid)=g{sky}.Z_SPEC(g{sky}.CATA_INDEX==cat{sky}.BCGRef(gid));
%         end
%     end
% end
% end
    

disp('loading shear...');
% e{1}=loadShearTBL(1,[1,2,40,41,42,53,57,58]);% ra,dec,chi1,chi2,sigma_e,z,z-sig,z+sig
% e{2}=loadShearTBL(2,[1,2,40,41,42,53,57,58]);
% e{3}=loadShearTBL(3,[1,2,40,41,42,53,57,58]);
% save shearmap.mat e

load shearmap
for sky=1:3
    e{sky}(:,5)=e{sky}(:,5)*2;
end  %change the error to error in chi (or e) rather than in epsilon


if isfield(macros,'ZGAMA')
if abs(macros.ZGAMA)>=2  %sdss photoz
    load src_sdss_match.mat nmatch idmatch zmatch zerrmatch
    if abs(macros.ZGAMA)==2  %sdss zcc2
        for i=1:3
            tmp=zmatch{i}(:,2)>-100;
            idtmp=find(nmatch{i}>0);
            e{i}(idtmp(tmp),6:8)=[zmatch{i}(tmp,2),zmatch{i}(tmp,2)-zerrmatch{i}(tmp,2),zmatch{i}(tmp,2)+zerrmatch{i}(tmp,2)];
        end
    end
     if abs(macros.ZGAMA)==3  %sdss zd1
         for i=1:3
             tmp=zmatch{i}(:,3)>-100;
             idtmp=find(nmatch{i}>0);
             e{i}(idtmp(tmp),6:8)=[zmatch{i}(tmp,3),zmatch{i}(tmp,3)-zerrmatch{i}(tmp,3),zmatch{i}(tmp,3)+zerrmatch{i}(tmp,3)];
         end
     end
     if abs(macros.ZGAMA)==4  %sdss template
         for i=1:3
             tmp=zmatch{i}(:,1)>-100;
             idtmp=find(nmatch{i}>0);
             e{i}(idtmp(tmp),6:8)=[zmatch{i}(tmp,1),zmatch{i}(tmp,1)-zerrmatch{i}(tmp,1),zmatch{i}(tmp,1)+zerrmatch{i}(tmp,1)];
         end
    end
end
if macros.ZGAMA>0
load GAMA_SRC_matchZ
for i=1:3
    e{i}(source_idmatch{i},6:8)=[gama_zmatch{i}(:,1),gama_zmatch{i}(:,1)-gama_zmatch{i}(:,2)/3e5,gama_zmatch{i}(:,1)+gama_zmatch{i}(:,2)/3e5];
end
end
end

cd /work/Projects/Lensing/output
%% photoz errors
% figure;
% plot(e{1}(:,6),e{1}(:,7),'.');hold on;
% plot(e{1}(:,6),e{1}(:,8),'r.');
% f=e{1}(:,6)>0;
% figure;linhist((e{1}(f,6)-e{1}(f,7))./(1+e{1}(f,6)),0:0.02:1,'stairs');  % cc2 peaks at sigma(z)=0.07*(1+z)
% figure;linhist((e{1}(f,8)-e{1}(f,6))./(1+e{1}(f,6)),0:0.02:1,'stairs');
% err=sqrt(((e{1}(f,6)-e{1}(f,7)).^2+(e{1}(f,8)-e{1}(f,6)).^2)/2)./(1+e{1}(f,6));
% figure;[x,y]=linhist(err,0:0.01:0.5,'stairs');  % cc2 peaks at sigma(z)=0.07*(1+z)
% title(num2str(x(y==max(y(x>0.02)))));
%% params
global OmegaM OmegaL 

OmegaM=0.25;OmegaL=0.75;

nbin=60; %radial bins
nbinm=4;nbinz=2;
%% linklist
global ll
ll=cell(3,1);
ngrids=[40,20];
if exist(['GAMAlinklist_',num2str(ngrids(1)),'_',num2str(ngrids(2)),'.mat'],'file')
    load(['GAMAlinklist_',num2str(ngrids(1)),'_',num2str(ngrids(2)),'.mat']);
else
    for sky=1:3
        ll{sky}.ngrid=ngrids;
        [ll{sky}.grids,ll{sky}.xrange,ll{sky}.yrange,ll{sky}.step]=linklist(e{sky}(:,1),e{sky}(:,2),ll{sky}.ngrid);
    end
    eval(['save GAMAlinklist_',num2str(ngrids(1)),'_',num2str(ngrids(2)),'.mat ll'])
end
for sky=1:3
    grp{sky}.HaloMass=sqrt(grp{sky}.LumMass.*grp{sky}.MassProxy);
end
%% initialize
labels=mark_mz_bin(macros.PROXY);
mr=zeros(nbin,nbinm,nbinz);   % n radial bins, nbinm mass bins, nbinz redshift bins
nr=mr;sr=mr;wr=mr;emr=mr;esr=sr;rr=mr;rr2=rr;sr_pred=sr;
ngrps_used=zeros(nbinm,nbinz);

rmin=macros.RMIN;
Mav=zeros(nbinm,nbinz,3);
for sky=1:3
for i=1:nbinm
    for j=1:nbinz
Mav(i,j,sky)=mean(grp{sky}.HaloMass(grp{sky}.mid==i&grp{sky}.zid==j))/1e10;
    end
end
end
Mav(isnan(Mav))=0;
if macros.RMAX<0
rmax=comoving_200b_radius(sum(Mav,3)./sum(logical(Mav),3),OmegaM)/1e3*(-1*macros.RMAX); %3*rv, Mpc/h
rmax=repmat(max(rmax,[],2),1,nbinz);
else
rmax=zeros(nbinm,nbinz)+macros.RMAX; % 10Mpc for all
end