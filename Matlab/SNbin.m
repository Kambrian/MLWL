function sn=SNbin(zl,M,N)
% M: in units 10^10Msun/h
% N: number of lenses at zl with mass M

n0=1.2; %source number density, arcmin^-2;
nstack=N*n0; %stacked number density
se=0.4; % square noise of ellipticity, including intrinsic and measurement, single component value
rmin=0.2;rmax=2;

[sigav,r200]=nfw_avDeltSig([rmin,rmax],M,zl);
theta200=r200/comoving_dist(0.3,0.7,zl)*(1+zl)/pi*180*60; %angular seperation in arcminutes
sn=sigav*sqrt(pi*(rmax^2-rmin^2)*theta200^2)*sqrt(nstack)*2/sigma_crit_eff(zl)/(se);

