function DL=lum_dist(OmegaM,z)
%luminosity dist fitting formula for flat cosmology
% in  units of Mpc/h
% accurate to within 0.3% for OmegaM=(0.2~1) and z=(0.01~inf).
%ref: M. Adachi and M. Kasai (arXiv:1012.2670)


x0=(1-OmegaM)/OmegaM;
F0=(2 + 2.641.*x0 + 0.8830.*x0.^2 + 0.05313.*x0.^3)./(1 + 1.392.*x0 + 0.5121.*x0.^2 + 0.03944.*x0.^3);
x=(1-OmegaM)/OmegaM./(1+z).^3;
Fx=(1+z).^(-1/2).*(2 + 2.641.*x + 0.8830.*x.^2 + 0.05313.*x.^3)./(1 + 1.392.*x + 0.5121.*x.^2 + 0.03944.*x.^3);
DL=3000*(1+z)./sqrt(OmegaM).*(F0-Fx);