function n=num_dens_prof_toy(r,z,zc,Mc)
%number density profile of galaxies around groups
%toy model: isothermal+uniform background+photoz blur+depth decay
%input: r, projected seperation (Mpc/h)
%           z, satellite redshift;
%               r and z must be the same size, forming pairs
%           zc, group redshift
%           Mc, group virial mass (10^10Msun/h)
%output: n, number density of galaxies at (r,z), (Mpc/h)^-3
%$

  
    G=43.0071; % gravity const
    H=100; %hubble param, km/s/(Mpc/h)
    virialF_b=200; %virial factor, 200 times background
    OmegaM=0.3;
    sigmz=0.1;

    rv=(Mc*2*G/virialF_b/OmegaM/H^2)^(1/3);
    req=sqrt(virialF_b/3)*rv;
%     rc=comoving_dist(OmegaM,1-OmegaM,zc);


n=num_dens_prof_blur(r,z).*exp(-8.8*z);
% n=r;
% for j=1:numel(r)
% n(j)=num_dens_prof_iso(r(j),z(j));
% end
    
 function Hz=Hubblez(z)
        Hz=H * sqrt(OmegaM .*(1+z).^3+(1-OmegaM));
 end

    function n=num_dens_prof_blur(r,z)
        n=z;
        for i=1:numel(z)
        n(i)=quadgk(@(zp) num_dens_prof_iso(r(i),zp).*exp(-(zp-z(i)).^2./2/sigmz^2)/sqrt(2*pi)/sigmz./Hubblez(zp),0,zc)+...
            quadgk(@(zp) num_dens_prof_iso(r(i),zp).*exp(-(zp-z(i)).^2./2/sigmz^2)/sqrt(2*pi)/sigmz./Hubblez(zp),zc,1.5);
        n(i)=n(i)*Hubblez(z(i));
        end
    end

    function n=num_dens_prof_iso(ri,z)
    %isothermal+background, without blur or decay
  
%     rz=comoving_dist(OmegaM,1-OmegaM,z)-rc;
   rz=AD_dist_flat(OmegaM,zc,z).*(1+z);
    
    d=sqrt(ri^2+rz.^2)/req;
    n(d<1)=d(d<1).^(-2);
%     n(d<1e-2)=1e4; % to avoid singularity at center
    n(d>=1)=1;
    end
end