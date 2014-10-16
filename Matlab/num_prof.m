function [data,datamock]=num_prof()
% to check the number density profile of the shear calc
clear global
% global macros
macros.PROXY='halomass';  %mass proxy for stacking
macros.LENSCATID=0; % catalogue used as lens: 0 means real GAMA, 1 means mockGAMA, <0 means random catalogue
macros.STACK_COMOV=1;  %comoving stacking
macros.RADIAL_LOGBIN=1;  %radial bin in logscale
macros.BCG_CENTER=1;     %use BCG as center
macros.MULT_MIN=3;      %minimum multiplicity of groups
macros.PHOTOZ_DENS_ERR2_MAX=0;  %maximum photoz uncertainty in SigmaCrit, set to 0 means do not calc nor include this error
macros.PHOTOZ_DENS_ERR_INC=0; %include photoz error in calc; PHOTOZ_DENS_ERR2_MAX must be >0 for this to take effect
macros.RMIN=0.02;  %rmin for radial bin, in Mpc/h
macros.RMAX=10;   %rmax for radial bin, in Mpc/h; if <0, set rmax=-RMAX*rvir;
macros.ZGAMA=1;
% macros.ZMAX;  %zmax of groups
if macros.PHOTOZ_DENS_ERR2_MAX==0
    macros.PHOTOZ_DENS_ERR_INC=0;  %disable zerr
end

ltyp={'ro-','ro--';'g>-','g>--';'bs-','bs--';'kd-','kd--';};
ldtyp={'ro:';'g>:';'bs:';'kd:';};
lltyp={'ro-';'g>-';'bs-';'kd-';};
mtyp={'ro';'g>';'bs';'kd';};
color={'r';'g';'b';'k'};
dir='/home/kam/Projects/Lensing/output/';

data=load_nprof(macros);
macromock=macros;
macromock.LENSCATID=-1;
% macromock.ZGAMA=1;
% macromock.ZDIFF=0.01;
datamock=load_nprof(macromock);

myfigure;
for i=1:4
    y=datamock.y{i};
    ey=datamock.ey{i};
%     y=y/nz_factor(datamock.zmean(i))*(pi/180/60)^2;
%     ey=ey/nz_factor(datamock.zmean(i))*(pi/180/60)^2;
    h1=errorbar(datamock.rvar{i}',y,ey,ldtyp{i});hold on;
    y=data.y{i};
    ey=data.ey{i};
%     y=y/nz_factor(data.zmean(i))*(pi/180/60)^2;
%     ey=ey/nz_factor(data.zmean(i))*(pi/180/60)^2;
    h2=errorbar(data.rvar{i}',y,ey,lltyp{i},'markerfacecolor',color{i});hold on;
end

set(gca,'xscale','log');
legend([h1,h2],'mock','gama');
xlabel('r/Mpc/h');
% ylabel('$<n>/(arcmin)^{-2}$');
% print('-depsc',[dir,'n_prof_cmp.eps']);
ylabel('$<n>/(Mpc/h)^{-2}$');set(gca,'yscale','log');
% print('-depsc',[dir,'n_prof_cmp_comv.eps']);
% print('-depsc','n_prof_mock.eps');
end

function data=load_nprof(macro)

global cat
nbinvar=[3,6,12,15];

switch macro.LENSCATID
    case 1
        cat=cell(1,3);
        catall=loadGAMAcsv('/home/kam/Projects/Lensing/data/mocksv2/mockgroupcatv2opt194lim194vol1.csv');
        [cat{1},cat{2},cat{3}]=split_mockcat(catall);
    otherwise
        cat{1}=loadGAMAcsv('/home/kam/Projects/Lensing/data/groupsv2/groupcatv2G09opt194lim194.csv');
        cat{2}=loadGAMAcsv('/home/kam/Projects/Lensing/data/groupsv2/groupcatv2G12opt194lim194.csv');
        cat{3}=loadGAMAcsv('/home/kam/Projects/Lensing/data/groupsv2/groupcatv2G15opt194lim194.csv');
%         error('unknown lens catalogue ');
end

for sky=1:3
    cat{sky}.HaloMass=sqrt(cat{sky}.LumMass.*cat{sky}.Mass);
end
labels=mark_mz_bin(macro.PROXY);
catall=structcat([cat{1},cat{2},cat{3}]);

data=gama_rebin(macro,nbinvar);
load(['/home/kam/Projects/Lensing/output/',gamaWL_fname(macro),'.mat']);
rmin=macro.RMIN;
data.rmax=rmax;
data.zmean=zeros(4,1);data.dl=zeros(4,1);data.v=cell(4,1);
f=catall.Mult>2;
for i=1:4
    data.zmean(i)=sqrt(mean(catall.MedianZ(f&catall.mid==i).^-2))\1;  % because nz_factor~z^-2, so the equivalent average is <z^-2>^-0.5;
    data.dl(i)=AD_dist_flat(0.3,0,data.zmean(i));
    tmp=logspace(log10(rmin),log10(rmax(i,1)),nbinvar(i)+1);
    data.v{i}=pi*diff(tmp.^2);
end

data.nmean=zeros(1,4);
data.y=cell(4,1);data.ey=cell(4,1);
for i=1:4
    data.nmean(i)=sum(data.nvar{i})/rmax(i,1)^2/pi/data.ngrps(i);
    data.nmean(i)=data.nmean(i)/nz_factor(data.zmean(i))*(pi/180/60)^2;
    data.y{i}=(data.nvar{i}'./data.v{i})/data.ngrps(i);
    data.ey{i}=sqrt(data.nvar{i}')./data.v{i}/data.ngrps(i);
end
end