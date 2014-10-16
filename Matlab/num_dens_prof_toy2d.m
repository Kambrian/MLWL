function n=num_dens_prof_toy2d(r,zc,Mc,dzmin,sigmz)
%projected number density profile of galaxies around groups
%toy model: isothermal+uniform background+photoz blur+depth decay
%input: r, projected seperation (Mpc/h)
%           zc, group redshift
%           Mc, group virial mass (10^10Msun/h)
%           dzmin, minimum redshift seperation to start the projection
%output: n, number density of galaxies at (r,z), (Mpc/h)^-3
%$

    c=3e5;
    G=43.0071; % gravity const
    H=100; %hubble param, km/s/(Mpc/h)
    virialF_b=200; %virial factor, 200 times background
    OmegaM=0.3;
%     sigmz=0.1;

    rv=(Mc*2*G/virialF_b/OmegaM/H^2)^(1/3);
    req=sqrt(virialF_b/3)*rv
    
    n=r;
    nb=quad(@(z) num_dens_prof_blur(10*req,z).*exp(-8.8*z)*c./Hubblez(z),zc+dzmin,1.5);

for i=1:numel(r)
n(i)=quad(@(z) num_dens_prof_blur(r(i),z).*exp(-8.8*z)*c./Hubblez(z),zc+dzmin,1.5);
end
n=n/nb;

    function Hz=Hubblez(z)
        Hz=H * sqrt(OmegaM .*(1+z).^3+(1-OmegaM));
    end

    function n=num_dens_prof_blur(ri,z)
        n=z;
        for j=1:numel(z)
%         n(j)=quadgk(@(zp) num_dens_prof_iso(ri,zp).*exp(-(z(j)-zp).^2./2/sigmz^2)/sqrt(2*pi)/sigmz./Hubblez(zp),0,1.5)*Hubblez(z(j));
        n(j)=quadgk(@(zp) num_dens_prof_iso(ri,zp).*exp(-(z(j)-zp).^2./2/sigmz^2)/sqrt(2*pi)/sigmz./Hubblez(zp),0,zc)+...
        quadgk(@(zp) num_dens_prof_iso(ri,zp).*exp(-(z(j)-zp).^2./2/sigmz^2)/sqrt(2*pi)/sigmz./Hubblez(zp),zc,1.5);
        n(j)=n(j)*Hubblez(z(j));
        end
    end

    function n=num_dens_prof_iso(ri,z)
    %isothermal+background, without blur or decay
  
    rz=AD_dist_flat(OmegaM,zc,z).*(1+z);
    
    d=sqrt(ri^2+rz.^2)/req;
    n(d<1)=d(d<1).^(-2);
%     n(d<1e-2)=1e4;
    n(d>=1)=1;
    end
end