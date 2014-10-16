tmp=importdata('/work/Projects/Lensing/data/srccat-gama/CFHTLens/CFHTLens.tsv');
data=tmp.data;
%%
ra=data(:,1);
dec=data(:,2);
area=diff(minmax(ra'))*diff(minmax(dec'))
%%
global cfht_z
filter=data(:,6)>0.1;
cfht_z=data(filter,6);
n0=numel(cfht_z)/area/3600; %3.09
%%
M=1e2;zl=1.3;
%n0=1.2; %source number density, arcmin^-2;
nlens=10000; %FASTSOUND Galaxies
nstack=nlens*n0; %stacked number density
errchi=0.4*sqrt(2); % noise of chi ellipticity, including intrinsic and measurement, two component value
rmin=0.2;rmax=2;

[sigav,r200]=nfw_avDeltSig([rmin,rmax],M,zl);
dl=comoving_function(comoving_dist(0.3,0.7,zl),1)/(1+zl);
r=logspace(log10(rmin*r200),log10(rmax*r200),6);  % logscale radial bin
% has good over S/N trend
% r=linspace((rmin*r200),(rmax*r200),20); % linear scale radial has highest inner S/N
sig=nfw_DeltSig(r,M,zl,0);
theta=r/dl/pi*180*60; %angular seperation in arcminutes
sn=sig(1:end-1).*sqrt(pi*diff(theta.^2))*sqrt(nstack)*2*sqrt(2)/sigma_crit_eff_cfht(zl)/errchi;
noise=sig(1:end-1)./sn;

[snav,snav_n,snmass]=SNmass_cfht(zl,M,nlens);
myfigure;semilogx(r(2:end),sn);hold on; 
plot(r,repmat(snav,size(r)),'r-');
plot(r,repmat(snav_n,size(r)),'g-');
plot(r,repmat(snmass,size(r)),'k-');
l=legend('$\Delta\Sigma_{wght}(r)$','$<\Delta\Sigma_{wght}>$','$<\Delta\Sigma>$','$M_{\zeta}$');
set(l,'interpreter','latex','location','west','box','off');
xlabel('R/(Mpc/h)');ylabel('S/N');
% print('-depsc','SN_FASTSOUND_CFHT.eps');
%%
figure;
errorbar(r(2:end)/r200,normrnd(sig(1:end-1),noise),noise,'o');hold('on');
plot(r(2:end)/r200,sig(1:end-1),'r-');
xlabel('r/r200');
ylabel('DSig');
set(gca,'xscale','log','yscale','log');
print('-depsc','MockSig4_FASTSOUND_CFHT.eps');