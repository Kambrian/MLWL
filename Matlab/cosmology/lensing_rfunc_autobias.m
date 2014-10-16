function [y,ynfw,ylin,rv]=lensing_rfunc_autobias(r,M,z,c)
% NFW+linear corr lensing profile
% r: comoving, Mpc/h
% M: 10^10Msun/h
% b: bias factor
% z: redshift. suggest to fix this as average z
% c: concentration. specify c<0 to automatically use theoretical M-c relation
% cosmology is fixed

OmegaM=0.3;H0=100;G=43.0071;
rhob=OmegaM*3*H0^2/8/pi/G;
virtype=2;
[ynfw,rv]=nfw_DeltSig(r./(1+z),M,z,virtype,c);
ynfw=ynfw/(1+z)^2;
rv=rv*(1+z);%convert to comoving
b=linear_bias(M,z,0.8);
ylin=b*lin_corr_DeltSig(r)*rhob*growth_factor(OmegaM,z)^2;
y=ynfw+ylin;
% figure;
% plot(r,y,'-');
% hold on;
% plot(r,ynfw,':');
% plot(r,ylin,'--');
% set(gca,'xscale','log','yscale','log');
% figure;
% plot(r,ylin./ynfw,'-');
% hold on;
% plot(r,ones(size(r)),'--');
% set(gca,'xscale','log','yscale','log');