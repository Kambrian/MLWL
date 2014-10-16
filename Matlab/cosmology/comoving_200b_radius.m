function R200b=comoving_200b_radius(M200b,OmegaM0)
% input: 
%       M200b: halo mass, in units of 10^10Msun/h
%       OmegaM0: z=0 matter density param
% output:
%       comoving virial radius when <rho(<r200b)>=200rho_m. 
%       in units of kpc/h
%       the same at any redshift.

G=43007.1;
R200b=(G*M200b/OmegaM0).^(1/3);