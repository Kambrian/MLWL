function ms=charact_mass(z,OmegaM,h,fb,sig8)
%return the characteristic mass ms in units 10^10M_sun/h
% sigm(ms)=delta_sc(z)
addpath('/work/Projects/HBT/code/v8.4/anal/Matlab/show');

w=collapse_barrier(1./(1+z),OmegaM);
ms=z;
for i=1:numel(z)
% ms(i)=fzero(@(m) mass_variance(m,OmegaM,h,fb,sig8)-w(i).^2,2e3);
logm=fzero(@(logm) mass_variance(exp(logm),OmegaM,h,fb,sig8)-w(i).^2,2);
ms(i)=exp(logm);
end