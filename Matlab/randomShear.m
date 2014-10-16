macro=default_WLparam();
macro.LENSCATID=-1; % catalogue used as lens: 0 means real GAMA, 1 means mockGAMA, <0 means random catalogue, >1 means random permutation catalogue
macro.RMAX=300;   %rmax for radial bin, in Mpc/h; if <0, set rmax=-RMAX*rvir;
macro.ZGAMA=2; %cc2
% macro.RTRIM=0; %trim radius
disp(macro)
Nrand=500;
%%
nbinvar=[30,30,30,30];
% nbinvar=[3,3,3,3];
nbinvar=[3,6,10,12];
randdata=cell(Nrand,1);
for i=1:Nrand
%     macro.LENSCATID=-i;
%     randdata{i}=gama_rebin(macro,nbinvar);
    file=['shearprof300_rand',num2str(i),'.dat'];
    randdata{i}=shearprof_rebin(file,nbinvar);
end

datar.nmean=cell(4,1);
datar.nstd=cell(4,1);
datar.sav=cell(4,1);
datar.smean=cell(4,1);  %average signal; you need 10^6 realizations to reduce smean from 10^3 to 1
datar.sstd=cell(4,1);
datar.esmean=cell(4,1);  %esmean is comparable to smean, means the uncertainty in the final estimation from the 50 measurements
datar.wmean=cell(4,1);
datar.wstd=cell(4,1);
for mbin=1:4
%mean and var in number profile
nvars=zeros(numel(randdata{1}.nvar{mbin}),Nrand);
for i=1:Nrand
nvars(:,i)=randdata{i}.nvar{mbin};
end
datar.nmean{mbin}=mean(nvars,2);
datar.nstd{mbin}=std(nvars,0,2);
%mean and var in shear signal
svars=nvars;wvars=nvars;esvars=nvars;
for i=1:Nrand
    svars(:,i)=randdata{i}.svar{mbin};
    wvars(:,i)=randdata{i}.wvar{mbin};
    esvars(:,i)=randdata{i}.esvar{mbin};
end
datar.smean{mbin}=sum(svars.*wvars,2)./sum(wvars,2);  %weighted average signal
datar.sav{mbin}=mean(svars,2);  %average signal, almost identical to smean
% datar.smean{mbin}=svars(:,1);  %weighted average signal
datar.sstd{mbin}=std(svars,0,2);  %scatter among random catalogues,or single measurement error
% datar.sstd{mbin}=sqrt(mean(svars.^2,2));  %scatter among random catalogues,or single measurement error
datar.esmean{mbin}=sqrt(sum((esvars.*wvars).^2,2))./sum(wvars,2);  %weighted average error est, or error on the estimation using 50 catlogues as a whole 
                        % this captures the uncertainty in shape, because
                        % of random tests of shape scatter; but does it
                        % capture photo-z uncertainty as well???? should be
                        % yes
datar.wmean{mbin}=mean(wvars,2);
datar.wstd{mbin}=std(wvars,0,2);
end

%real data
macro.LENSCATID=0;
% data=shearprof_rebin('shearprof_sky0_quarter.dat',nbinvar);
data=shearprof_rebin('shearprof300.dat',nbinvar);
% data=gama_rebin(macro,nbinvar);
% load(['/home/kam/Projects/Lensing/output/',gamaWL_fname(macro),'.mat']);
rmin=macro.RMIN;
rmax=macro.RMAX;
data.rmax=macro.RMAX;
data.v=cell(4,1);
for i=1:4
    tmp=logspace(log10(rmin),log10(rmax),nbinvar(i)+1);
    data.v{i}=pi*diff(tmp.^2);
end
%%
data.y=cell(4,1);data.ey=cell(4,1);
datar.y=cell(4,1);datar.ey=cell(4,1);datar.eyp=cell(4,1);
data.wy=cell(4,1);
data.ewy=cell(4,1);
datar.wy=cell(4,1);datar.ewy=cell(4,1);datar.ewp=cell(4,1);
for i=1:4
    data.y{i}=(data.nvar{i}./data.v{i})/data.ngrps(i);
    data.ey{i}=sqrt(data.nvar{i})./data.v{i}/data.ngrps(i);
    data.wy{i}=(data.wvar{i}./data.v{i})/data.ngrps(i);
    data.ewy{i}=data.ey{i}.*data.wy{i}./data.y{i};
    datar.y{i}=(datar.nmean{i}'./data.v{i})/data.ngrps(i);
    datar.ey{i}=datar.nstd{i}'./data.v{i}/data.ngrps(i); %real scatter, gives the error for each measurement
    datar.eyp{i}=sqrt(datar.nmean{i}'*Nrand)/Nrand./data.v{i}/data.ngrps(i);%total poisson error for the 50 measurements, or error on the estimated mean
    datar.wy{i}=(datar.wmean{i}'./data.v{i})/data.ngrps(i);
    datar.ewp{i}=datar.eyp{i}.*datar.wy{i}./datar.y{i};
    datar.ewy{i}=(datar.wstd{i}'./data.v{i})/data.ngrps(i);
end
%%
ltyp={'ro-','ro--';'g>-','g>--';'bs-','bs--';'kd-','kd--';};
ldtyp={'ro:';'g>:';'bs:';'kd:';};
lltyp={'ro-';'g>-';'bs-';'kd-';};
mtyp={'ro';'g>';'bs';'kd';};
color={'r';'g';'b';'k'};
dir='/home/kam/Projects/Lensing/output/';
%% number profiles
myfigure;
for i=1:4
    y=datar.y{i};
    ey=datar.ey{i};
    h1=errorbar(data.rvar{i}',y,ey,ldtyp{i},'color','r');hold on;
    y=datar.y{i};
    ey=datar.eyp{i};
    h2=errorbar(data.rvar{i}',y,ey,mtyp{i},'color','g');hold on;
    y=data.y{i};
    ey=data.ey{i};
    h3=errorbar(data.rvar{i}',y,ey,lltyp{i},'color','b','markerfacecolor','b');hold on;
end
set(gca,'xscale','log');
l=legend([h1,h2,h3],'random+single err','rand+poisson','gama');
set(l,'location','southwest');
xlabel('r/(Mpc/h)');
ylabel('$<n>/(Mpc/h)^{-2}$');set(gca,'yscale','log');
% print('-depsc',[dir,'n_prof_cmp_rand.eps']);
% print('-depsc',[dir,'ntrim_prof_cmp_rand_300.eps']);

myfigure;
for i=1:4
    h1=errorbar(data.rvar{i}',datar.wy{i},datar.ewy{i},ldtyp{i},'color','g');hold on;
    h2=errorbar(data.rvar{i}',datar.wy{i},datar.ewp{i},mtyp{i},'color','b','markerfacecolor','b');hold on;
    h3=errorbar(data.rvar{i}',data.wy{i},data.ewy{i},lltyp{i},'color','r');hold on;
end
set(gca,'xscale','log');
l=legend([h1,h2,h3],'random+single err','rand+poisson','gama');
set(l,'location','southwest');
xlabel('r/(Mpc/h)');
ylabel('$<w>/(Mpc/h)^{-2}$');set(gca,'yscale','log');
% print('-depsc',[dir,'n_prof_cmp_rand.eps']);
% print('-depsc',[dir,'wtrim_prof_cmp_rand_300.eps']);

myfigure;
for i=1:4
    eyy=sqrt((datar.eyp{i}./datar.y{i}).^2+(data.ey{i}./data.y{i}).^2).*data.y{i}./datar.y{i};
    errorbar(data.rvar{i},data.y{i}./datar.y{i},eyy,lltyp{i});
    hold on;
    errorbar(data.rvar{i},datar.y{i}./mean(datar.y{i}),datar.eyp{i}./mean(datar.y{i}),ldtyp{i});
%     plot(data.rvar{i},data.wy{i}./datar.wy{i},mtyp{i});
end
plot([1e-2,100],[1,1],'--');
ylim([0,5]);
set(gca,'xscale','log');
xlabel('r/(Mpc/h)');ylabel('$n_{GAMA}/n_{rand}$');
% print('-depsc',[dir,'n_prof_rat_rand.eps']);
% print('-depsc',[dir,'ntrim_prof_rat_rand_300.eps']);

myfigure;
for i=1:4
    eyy=sqrt((datar.ewp{i}./datar.wy{i}).^2+(data.ewy{i}./data.wy{i}).^2).*data.wy{i}./datar.wy{i};
    errorbar(data.rvar{i},data.wvar{i}./datar.wmean{i}',eyy,lltyp{i});
    hold on;
%     errorbar(data.rvar{i},datar.y{i}./mean(datar.y{i}),datar.eyp{i}./mean(datar.y{i}),ldtyp{i});
end
plot([1e-2,100],[1,1],'--');
ylim([0,5]);
set(gca,'xscale','log');
xlabel('r/(Mpc/h)');ylabel('$w_{GAMA}/w_{rand}$');
% print('-depsc',[dir,'wtrim_prof_rat_rand_300.eps']);
%% error comparison
myfigure;
for l=1:4
% h1=plot(data.rvar{l},signlog(data.svar{l}-datar.smean{l}'),lltyp{l},'color','r','markerfacecolor','r');hold on;
h2=plot(data.rvar{l},(data.esvar{l}),lltyp{l},'markerfacecolor',color{l});hold on;
% h3=plot(data.rvar{l},signlog((datar.smean{l})),ldtyp{l},'color','b');hold on;
% h4=plot(data.rvar{l},signlog((datar.sav{l})),mtyp{l},'color','g');hold on;
h4=plot(data.rvar{l},(datar.sstd{l}),mtyp{l},'markersize',8);hold on;
% h5=plot(data.rvar{l},signlog(datar.esmean{l}),mtyp{l},'color','k');hold on;
end
% plot([1e-2,1e1],10.^[4.5,1.5],'c:');
% plot([1e-2,1e1],[-4.5,-1.5],'c:',[1e-2,1e1],[4.5,1.5],'c:');
legend([h2,h4],'gama err','rand scatter');
% legend([h1,h2,h3,h4,h5],'gama','gama err','rand','rand scatter','rand
% err');
% legend([h2,h3,h4,h5],'gama err','rand','rand scatter','randerr');
% legend([h1,h3],'rand','randerr');
xlabel('r/(Mpc/h)');ylabel('$\sigma_{\Delta\Sigma}$');
set(gca,'xscale','log');
set(gca,'yscale','log');
% print('-depsc',[dir,'wlErr_rand_err.eps']);
%%
figure;
handles=zeros(4,1);
for i=1:4
%     s=(data.svar{i}-datar.smean{i}').*data.wy{i}./datar.wy{i};
    s=(data.svar{i}-datar.smean{i}').*data.wy{i}./datar.wy{i};
    es=sqrt((data.esvar{i}.^2+datar.esmean{i}'.^2)./s.^2+data.ewy{i}.^2./data.wy{i}.^2+datar.ewy{i}.^2./datar.wy{i}.^2).*abs(s);
%     h=errorbar(data.rvar{i},signlog(s),signlog(s)-signlog(s-es),signlog(s+es)-signlog(s),mtyp{i});
   h=plot(data.rvar{i},signlog(s),mtyp{i});
   hold on;
   handles(i)=h(1);
end
plot([1e-2,1e1],[-4.5,-1.5],'c:',[1e-2,1e1],[4.5,1.5],'c:');

GAMA_WL_init;
catall=structcat([cat{1},cat{2},cat{3}]);
f=catall.Mult>=macro.MULT_MIN;
zmean=zeros(4,1);halomass_mean=zmean;halomass_err=zmean;
for i=1:4
zmean(i)=mean(catall.MedianZ(f&catall.mid==i));
halomass_mean(i)=mean(catall.HaloMass(f&catall.mid==i));
halomass_err(i)=std(catall.HaloMass(f&catall.mid==i));
end
M=halomass_mean/1e10;eM=halomass_err/1e10;
z1=zmean(1);
z2=zmean(2);      
z3=zmean(3);     
z4=zmean(4);  
r1=AD_dist_flat(0.3,0,z1)*1/180*pi*(1+z1);
r2=AD_dist_flat(0.3,0,z2)*1/180*pi*(1+z2);
r3=AD_dist_flat(0.3,0,z3)*1/180*pi*(1+z3);
r4=AD_dist_flat(0.3,0,z4)*1/180*pi*(1+z4);
plot([r1,r1],[-5,5],'color',color{1});
plot([r2,r2],[-5,5],'color',color{2});
plot([r3,r3],[-5,5],'color',color{3});
plot([r4,r4],[-5,5],'color',color{4});

% comoving density in comoving radius
for i=1:numel(M)
    z=zmean(i);
    plot(data.rvar{i},signlog(nfw_DeltSig(data.rvar{i}./(1+z),M(i),z,2)/(1+z)^2),'-','color',color{i});
    hold on;
end
% set(gca,'yscale','log');
set(gca,'xscale','log');
xlabel('r/(Mpc/h) comov');ylabel('slog(|\Delta\Sigma(r)|)');    %/(10^{10}M_\odot/h/(Mpc/h)^2)
% title('residual shear subtracted, dilution corrected');
title('trimmed gama+trimmed src')
% title('sky2');
% print('-depsc',[dir,'wlSig_300_doubletrim.eps']);
% print('-depsc',[dir,'wlSig_300_subtraction_amplication.eps']);
%% S/N
figure;
handles=zeros(4,1);
for i=1:4
%     s=(data.svar{i}-datar.smean{i}').*data.y{i}./datar.y{i};
    s=abs(data.svar{i})./data.esvar{i};
    h=plot(data.rvar{i},(s),mtyp{i});
%    h=plot(data.rvar{i},signlog(s),mtyp{i});
   hold on;
   handles(i)=h(1);
end
set(gca,'xscale','log')
%%
myfigure;
handles=zeros(4,1);
for i=1:4
    serr=data.svar{i};
  h= plot(data.rvar{i},datar.sstd{i}./data.esvar{i}',mtyp{i});
   hold on;
   handles(i)=h(1);
end
set(gca,'xscale','log');
xlabel('r/(Mpc/h) comov');ylabel('$\sigma_{rand}/\sigma_{gama}$');    %/(10^{10}M_\odot/h/(Mpc/h)^2)
print('-depsc','rand_noise_rat.eps');
% title('boosted');

myfigure;
handles=zeros(4,1);
for i=1:4
    serr=data.svar{i};
  h= plot(data.rvar{i},datar.smean{i}./data.svar{i}',mtyp{i});
   hold on;
   handles(i)=h(1);
end
set(gca,'xscale','log');hold on;plot(get(gca,'xlim'),[1,1]);
xlabel('r/(Mpc/h) comov');ylabel('$\Delta\Sigma_{rand}/\Delta\Sigma_{gama}$');    %/(10^{10}M_\odot/h/(Mpc/h)^2)