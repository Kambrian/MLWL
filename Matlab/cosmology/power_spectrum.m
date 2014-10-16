function p=power_spectrum(OmegaM,h,fb,sig8,k,z)
% linear matter power spectrum at z
% OmegaM,h,fb,sig8: 
%      OmegaM, hubble param, OmegaBaryon/OmegaM, sigma_8; all at z=0.
% k: (Mpc/h)^-1, comoving
% z: redshift 
% p: (Mpc/h)^3

ns=1;  %primordial index

wx=@(x) 3*(sin(x)-x.*cos(x))./x.^3;

delt8=quadgk(@(k) wx(k*8).^2.*k.^(2+ns).*TF_BBKS(OmegaM,h,fb,k).^2,0,inf,'RelTol',1e-2);  %this integral is contributed primarily by k from 0.1 to 1
A=sig8.^2/delt8*2*pi^2;

p=A*k.^ns.*TF_BBKS(OmegaM,h,fb,k).^2;

p=p*growth_factor(OmegaM,z).^2;