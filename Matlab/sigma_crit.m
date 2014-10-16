function sig=sigma_crit(OmgM,OmgL,zl,zs)
%critical density for lensing

c=3e5; %light speed, km/s
G=43.0071; %gadget unit modified from kpc/h to Mpc/h
% Hubble=100; %km/s/(Mpc/h)

if OmgM+OmgL==1
DL=AD_dist_flat(OmgM,0,zl);
DS=AD_dist_flat(OmgM,0,zs);
DLS=AD_dist_flat(OmgM,zl,zs);
else
chiL=comoving_dist(OmgM,OmgL,zl);
chiS=comoving_dist(OmgM,OmgL,zs);
DL=comoving_function(chiL,OmgM+OmgL)./(1+zl);
DS=comoving_function(chiS,OmgM+OmgL);  % (1+zs) factor ommitted
DLS=comoving_function(chiS-chiL,OmgM+OmgL); %(1+zs) factor ommitted
end

sig=c^2/4/pi/G.*DS./DL./DLS;   %10^10Msun/h/(Mpc/h)^2
