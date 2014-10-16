function DA=AD_dist_flat(OmegaM,z1,z2)
% DA=AD_dist_flat(OmegaM,z1,z2)
%generalized angular dist fitting formula for flat cosmology
% DA=scaleF(z2)*comoving_dist(z1,z2)
% in  units of Mpc/h
% accurate to within 0.3% for OmegaM=(0.2~1) and z=(0.01~inf).
%ref: M. Adachi and M. Kasai (arXiv:1012.2670; 1111.6396)
% set z1=0 for usual DA


x1=(1-OmegaM)/OmegaM./(1+z1).^3;
Fx1=(1+z1).^(-1/2).*(2 + 2.641.*x1 + 0.8830.*x1.^2 + 0.05313.*x1.^3)./(1 + 1.392.*x1 + 0.5121.*x1.^2 + 0.03944.*x1.^3);
x2=(1-OmegaM)/OmegaM./(1+z2).^3;
Fx2=(1+z2).^(-1/2).*(2 + 2.641.*x2 + 0.8830.*x2.^2 + 0.05313.*x2.^3)./(1 + 1.392.*x2 + 0.5121.*x2.^2 + 0.03944.*x2.^3);
DA=3000./(1+z2)./sqrt(OmegaM).*(Fx1-Fx2);