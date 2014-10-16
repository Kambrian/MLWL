function y=lensing_rfunc(r,M,b,z,c)
% NFW+linear corr lensing profile (comoving DSig)
% r: comoving, Mpc/h
% M: 10^10Msun/h
% b: bias factor
% z: redshift. suggest to fix this as average z
% c: concentration. specify c<0 to automatically use theoretical M-c relation
% cosmology is fixed

OmegaM=0.3;H0=100;G=43.0071;
rhob=OmegaM*3*H0^2/8/pi/G;
virtype=2;
y=nfw_DeltSig(r./(1+z),M,z,virtype,c)/(1+z)^2+b*lin_corr_DeltSig(r)*rhob*growth_factor(OmegaM,z)^2;