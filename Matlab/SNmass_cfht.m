function [snav,snav_normal,snmass]=SNmass_cfht(zl,M,N)
% M: in units 10^10Msun/h
% N: number of lenses at zl with mass M

global cfht_z
area=22.7150;
n0=numel(cfht_z)/area/3600;
%n0=3.09; %source number density, arcmin^-2;
nstack=N*n0; %stacked number density
errchi=0.4*sqrt(2); % noise of chi ellipticity, including intrinsic and measurement, two component value
erreps=errchi/2; % noise of epsilon ellipticity, two component
rmin=0.2;rmax=2;

dl=comoving_function(comoving_dist(0.3,0.7,zl),1)/(1+zl);
sigmeff=sigma_crit_eff_cfht(zl);
sigmeff_normal=sigma_crit_eff_normal_cfht(zl);
[sigav,r200]=nfw_avDeltSig([rmin,rmax],M,zl);
theta200=r200/dl/pi*180*60; %angular seperation in arcminutes
snmass=M*sqrt(nstack)/dl/r200/erreps/sigmeff_normal*180*60/pi^(3/2)*sqrt(3/2); % signal-to-noise for M_{zeta} from rvir to 2rvir
snav=sigav*sqrt(pi*(rmax^2-rmin^2)*theta200^2)*sqrt(nstack)*2*sqrt(2)/sigmeff/errchi; % signal to noise for average shear from rmin to rmax, Madabaulm estimator
snav_normal=snav*sigmeff/sigmeff_normal; % average shear, simple estimator

