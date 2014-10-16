macro=default_WLparam();
macro.LENSCATID=0; % catalogue used as lens: 0 means real GAMA, 1 means mockGAMA, <0 means random catalogue, >1 means random permutation catalogue
macro.RMAX=10;   %rmax for radial bin, in Mpc/h; if <0, set rmax=-RMAX*rvir;
macro.PROXY='lummass';
macro.ZGAMA=0;
disp(macro)
GAMA_WL_init;
%%
file=gamaWL_fname(macro);
load('-mat',[file,'.mat']);
%%
catall=structcat([grp{1},grp{2},grp{3}]);
f=catall.Mult>=macro.MULT_MIN;
zmean=zeros(4,1);halomass_mean=zmean;halomass_err=zmean;
for i=1:4
zmean(i)=mean(catall.Zfof(f&catall.mid==i));
halomass_mean(i)=mean(catall.HaloMass(f&catall.mid==i));
halomass_err(i)=std(catall.HaloMass(f&catall.mid==i));
end
M=halomass_mean/1e10;eM=halomass_err/1e10;
%%
ltyp={'ro-','ro--';'g>-','g>--';'bs-','bs--';'kd-','kd--';};
ldtyp={'ro:';'g>:';'bs:';'kd:';};
lltyp={'ro-';'g>-';'bs-';'kd-';};
mtyp={'ro';'g>';'bs';'kd';};
color={'r';'g';'b';'k'};
%% M_\zeta
myfigure;
% for j=1:size(mr,3)
%     for i=1:size(mr,2)
%         semilogy(r(:,i,j),mr(:,i,j),ltyp{i,j});
%         hold on;
%     end
% end
ltyp={'ro';'g>';'bs';'kd';};
for i=1:size(m,2)
        errorbar(rl(:,i),m(:,i),em(:,i),ldtyp{i});
        hold on;
end
set(gca,'xscale','log','yscale','log');
xlabel('r/(Mpc/h) comov');ylabel('$M_\zeta(<r)/(10^{10}M_\odot/h)$');    
h=legend(labels{1},labels{2},labels{3},labels{4});
set(h,'location','southeast','interpreter','latex');
%comoving fitting in logspace
% M=[5.4e3,1e3,7e3,4.6e4];  
% eM=[2e3,3e2,2.6e3,1.3e4];
% %physical fitting in logspace
% M=[7e3,1.5e3,5.7e3,1.4e5];
% eM=[3e3,7e2,1.7e3,5e4];
for i=1:4
plot(rl(:,i),repmat(M(i),size(r(:,i))),'-','color',color{i});
hold on;
plot(rl(:,i),repmat(M(i)-eM(i),size(r(:,i))),':','color',color{i});
plot(rl(:,i),repmat(M(i)+eM(i),size(r(:,i))),':','color',color{i});
end

% print('-depsc',['Mwl_',bintype,'_comov.eps']);
%%
myfigure;
% for j=1:size(sr,3)
%     for i=1:size(sr,2)
%         semilogy(r(:,i,j),sr(:,i,j),ltyp{i,j});
%         hold on;
%     end
% end
for i=1:size(s,2)
%         errorbar(r(:,i),s(:,i),es(:,i),ltyp{i});
%         errorbar(r(:,i),s(:,i),min(es(:,i),s(:,i)-1),es(:,i),ltyp{i});
        semilogx(r(:,i),signlog(s(:,i)),mtyp{i});
        hold on;
end
plot([1e-2,1e1],[-4.5,-1.5],'k--',[1e-2,1e1],[4.5,1.5],'k--');
%comoving fitting in logspace
% M=[5.4e3,1e3,7e3,4.6e4];
% eM=[2e3,3e2,2.6e3,1.3e4];
% %physical fitting in logspace
% M=[7e3,1.5e3,5.7e3,1.4e5];
% eM=[3e3,7e2,1.7e3,5e4];
% r=0.1:0.5:3.5;
%comoving density in comoving radius
for i=1:numel(M)
    z=zmean(i);
    plot(r(:,i),log10(nfw_DeltSig(r(:,i)./(1+z),M(i),z,2)/(1+z)^2),'-','color',color{i});
    hold on;
end
% set(gca,'yscale','log');
set(gca,'xscale','log');
xlabel('r/(Mpc/h) comov');ylabel('$\Delta\Sigma(r)/(10^{10}M_\odot/h/(Mpc/h)^2)$');    
h=legend(labels{1},labels{2},labels{3},labels{4});
set(h,'location','southeast','interpreter','latex');
% print('-depsc',['Swl_',bintype,'_comov.eps']);
%%
cftool(r(:,4),s(:,4),es(:,4).^-2)
cftool(r(:,3),s(:,3),es(:,3).^-2)
cftool(r(:,2),s(:,2),es(:,2).^-2)
cftool(r(:,1),s(:,1),es(:,1).^-2)
%%
for i=1:4
f=s(:,i)>0;
cftool(r(f,i),log(s(f,i)),s(f,i).^2./es(f,i).^2)
end
%comoving fitting in logspace
% M=[5.4e3,1e3,7e3,4.6e4];
% eM=[2e3,3e2,2.6e3,1.3e4];
% %physical fitting in logspace
% M=[7e3,1.5e3,5.7e3,1.4e5];
% eM=[3e3,7e2,1.7e3,5e4];

%% rebin
% nbinvar=[30,30,30,30];
nbinvar=[3,6,12,15];
% nbinvar=[6,12,20,30];
data=gama_rebin(macros,nbinvar);
% data=a;

myfigure;

handles=zeros(4,1);
for i=1:4
%    h=ploterr(rvar{i}',svar{i}',sqrt(ervar{i}-rvar{i}.^2)',esvar{i}',ltyp{i},'logx');
%    h=semilogx(data.rvar{i}',signlog(data.svar{i}'),mtyp{i},'markerfacecolor',color{i}); hold on;
%    h=semilogx(data.rvar{i}',signlog(data.esvar{i}'),mtyp{i});     hold on;
   h=errorbar(data.rvar{i},signlog(data.spredvar{i}),signlog(data.spredvar{i})-signlog(data.spredvar{i}-data.esvar{i}),signlog(data.spredvar{i}+data.esvar{i})-signlog(data.spredvar{i}),mtyp{i});
   hold on;
   handles(i)=h(1);

%         errorbar(r(:,i),s(:,i),min(es(:,i),s(:,i)-1),es(:,i),ltyp{i});
end
plot([1e-2,1e1],[-4.5,-1.5],'c:',[1e-2,1e1],[4.5,1.5],'c:');
%comoving fitting in logspace
% M=[5.4e3,1e3,7e3,4.6e4];
% eM=[2e3,3e2,2.6e3,1.3e4];
% %physical fitting in logspace
% M=[7e3,1.5e3,5.7e3,1.4e5];
% eM=[3e3,7e2,1.7e3,5e4];
% r=0.1:0.5:3.5;
% z=0;
% comoving density in comoving radius
for i=1:numel(M)
    z=zmean(i);
    plot(data.rvar{i},log10(nfw_DeltSig(data.rvar{i}./(1+z),M(i),z,2)/(1+z)^2),'-','color',color{i});
    hold on;
end
% set(gca,'yscale','log');
set(gca,'xscale','log');
xlabel('r/(Mpc/h) comov');ylabel('$slog(|\Delta\Sigma(r)|)$');    %/(10^{10}M_\odot/h/(Mpc/h)^2)
h=legend(handles,labels{1},labels{2},labels{3},labels{4});
set(h,'location','southeast','interpreter','latex','box','off','color','none','fontsize',12);
% print('-depsc','wlSig_PhotozFilterdErrCons.eps');
% print('-depsc','wlSig_AllErrCons.eps');
% print('-depsc','wlSig_NewZMultErrCons.eps');
% print('-depsc','wlSig_BCGNewZMultErrCons.eps');
% print('-depsc','wlSigLum_BCGNewZMultErrCons.eps');
% print('-depsc','wlSigDyn_BCGNewZMultErrCons.eps');
% print('-depsc','wlSigBCGMag_BCGNewZMultErrCons.eps');
% print('-depsc','wlSigComvsize_BCGNewZMultErrCons.eps');
%%
nbinvar=[3,6,12,15];
macros.LENSCATID=0;
macros.MULT_MIN=2; 
macros.PHOTOZ_DENS_ERR2_MAX=0; %not filtering
macros.PHOTOZ_DENS_ERR_INC=0; %not counting z-err

macro=repmat(macros,[2,2,2]);
%z-filtering
for i=1:2
    for j=1:2
macro(1,i,j).PHOTOZ_DENS_ERR2_MAX=0;
macro(2,i,j).PHOTOZ_DENS_ERR2_MAX=2;
    end
end
%counting z-err
for i=1:2
    for j=1:2
macro(i,1,j).PHOTOZ_DENS_ERR_INC=0; 
macro(i,2,j).PHOTOZ_DENS_ERR_INC=1; 
    end
end
%Mult-filtering
for i=1:2
    for j=1:2
macro(i,j,1).MULT_MIN=2; 
macro(i,j,2).MULT_MIN=3;  
    end
end
macro(1,2,1).PHOTOZ_DENS_ERR2_MAX=100; %no z-filter but incl z-err
macro(1,2,2).PHOTOZ_DENS_ERR2_MAX=inf; %no z-filter but incl z-err

data=cell(2,2,2);
for i=1:2
    for j=1:2
        for k=1:2
            data{i,j,k}=gama_rebin(macro(i,j,k),nbinvar);
        end
    end
end

myfigure;
for l=1:4
    h1=plot(data{1,1,1}.rvar{l},signlog(data{1,1,1}.svar{l}),ldtyp{l});hold on;
    h2=plot(data{1,1,2}.rvar{l},signlog(data{1,1,2}.svar{l}),lltyp{l},'markerfacecolor',color{l});hold on;
end
h=legend([h1,h2],'all','$Mult>3$');set(h,'interpreter','latex','location','southeast');
xlabel('r/(Mpc/h) comov');ylabel('$sgn*log\{abs[\Delta\Sigma(r)]\}$');
plot([1e-2,1e1],[-4.5,-1.5],'c:',[1e-2,1e1],[4.5,1.5],'c:');
set(gca,'xscale','log');
print('-depsc','wlSig_cmp__MultFilter.eps');

%the filter in photoz can supress error estimation, but at the expense of
%messing the signal up a little, because your photo-z error estimation may
%be rubbish, so the filtering would be rubbish and in fact u do not get any
%better measurement from photo-z filtering
for i=1:2
    for j=1:2
figure;
for l=1:4
plot(data{1,i,j}.rvar{l},signlog(data{1,i,j}.esvar{l}),lltyp{l});hold on;
plot(data{2,i,j}.rvar{l},signlog(data{2,i,j}.esvar{l}),ldtyp{l},'markerfacecolor',color{l});hold on;
end
plot([1e-2,1e1],[-4.5,-1.5],'c:',[1e-2,1e1],[4.5,1.5],'c:');
set(gca,'xscale','log');
    end
end

%Mult>=3 gives better signal than Mult>=2, but of course increase
%statistical error slightly
for i=1
    for j=1:2
figure;
for l=1:4
plot(data{1,i,j}.rvar{l},signlog(data{1,i,j}.svar{l}),lltyp{l});hold on;
plot(data{2,i,j}.rvar{l},signlog(data{2,i,j}.svar{l}),ldtyp{l},'markerfacecolor',color{l});hold on;
end
plot([1e-2,1e1],[-4.5,-1.5],'c:',[1e-2,1e1],[4.5,1.5],'c:');
set(gca,'xscale','log');
    end
end


macromock=macros;
macromock.LENSCATID=1;
macromock.MULT_MIN=3; 
macromock.PHOTOZ_DENS_ERR2_MAX=2; %filtering
macromock.PHOTOZ_DENS_ERR_INC=1; %counting z-err
datamock=gama_rebin(macromock,nbinvar);
macromock2=macromock;
macromock2.PHOTOZ_DENS_ERR2_MAX=inf; %filtering
macromock3=macromock;
macromock3.PHOTOZ_DENS_ERR2_MAX=0; %filtering
macromock3.PHOTOZ_DENS_ERR_INC=0; %counting z-err
datamock2=gama_rebin(macromock2,nbinvar);
datamock3=gama_rebin(macromock3,nbinvar);
% the inclusion of photo-z error overestimates the fluctuation
% however, the shear-error along tend to over-estimate the error in the
% inner part for the largest mass bin
figure;
for l=1:4
plot(datamock.rvar{l},signlog(datamock.svar{l}),lltyp{l});hold on;
plot(datamock.rvar{l},signlog(datamock.esvar{l}),mtyp{l});hold on;
end
plot([1e-2,1e1],[-4.5,-1.5],'c:',[1e-2,1e1],[4.5,1.5],'c:');
set(gca,'xscale','log');
figure;
for l=1:4
plot(datamock2.rvar{l},signlog(datamock2.svar{l}),lltyp{l});hold on;
plot(datamock2.rvar{l},signlog(datamock2.esvar{l}),mtyp{l});hold on;
end
plot([1e-2,1e1],[-4.5,-1.5],'c:',[1e-2,1e1],[4.5,1.5],'c:');
set(gca,'xscale','log');
print('-depsc','wlErr_mock__all_ErrZinc.eps');
figure;
for l=1:4
plot(datamock3.rvar{l},signlog(datamock3.svar{l}),lltyp{l});hold on;
plot(datamock3.rvar{l},signlog(datamock3.esvar{l}),mtyp{l});hold on;
plot(data{1,1,2}.rvar{l},signlog(data{1,1,2}.esvar{l}),mtyp{l},'markerfacecolor',color{l});hold on;
end
plot([1e-2,1e1],[-4.5,-1.5],'c:',[1e-2,1e1],[4.5,1.5],'c:');
set(gca,'xscale','log');
print('-depsc','wlErr_mock__all_ErrCons.eps');
%    h=semilogx(data.rvar{i}',signlog(data.esvar{i}'),mtyp{i});     hold
%    on;
%% S/N
figure;
for i=1:4
   h=ploterr(data.rvar{i}',data.svar{i}'./data.esvar{i}',sqrt(data.ervar{i}-data.rvar{i}.^2)',0,lltyp{i},'logx');   hold on;
   handles(i)=h(1);
end
%% compare our bin abundance with that of Johnston et al. 2007
myfigure;
[x,y,dy]=loghist(catall.HaloMass,12,'stairs','k');hold on;
[x,y,dy]=loghist(catall.HaloMass(catall.Mult>2),12,'stairs','k--');
xlabel('Group Mass');
ylabel('dN');
hold on;
x=[4.3,5.3,8,13,9.7,12.7,25.5,42.3,74.5,123,199,503]*1e12;
y=[58788,27083,14925,8744,5630,3858,6196,4427,1711,787,272,47];
% ex=[0.45,0.65,1.34,1.65,2.28,3.36,2.86,3.42,7.46,11.28,24.81,87.61]*1e12;
hold on;plot(x,y,'mo--');
x=halomass_mean;
y=sum(ngrps_used,2);
hold on;plot(x,y,'rs--');
legend('GAMA','Mult>2','Jhonston','MassBin','location','best');
print('-depsc','GroupBinCmp.eps');