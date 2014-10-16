function chi=comoving_dist(OmgM,OmgL,z)
% comoving distance coordinate 
% output:
%    chi=int(dr/sqrt(1-k*r^2),0,r)----unit:Mpc/h, with H0=100km/s/Mpc*h
% input: 
%    OmgM, OmgL-----cosmological params at z=0
%           z  -----target redshift 


dDdz=@(z) 3000./sqrt((1+z).^2.*(1+OmgM*z)-z.*(2+z)*OmgL);  %c/H0=3000Mpc/h

chi=z;
for i=1:numel(z)
chi(i)=quad(dDdz,0,z(i));
end