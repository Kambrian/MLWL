function [sig,rv]=sub_surf_overdensity(r,M,z,virtype,cn)
%Sigma=int(rho_sub,z,-inf,inf)/rho_background for sub halo number density(Physical),
% surface number overdensity of subhalos for a NFW halo sunperimposed onto a uniform background
% suppose subhalo number density follows dark matter density from Rv and
% beyond, and follows einasto profile with alpha=1 within Rv.
% rv: physical virial radius
%
% r: Mpc/h, physical
% M: 10^10Msun/h
% z: redshift
% virtype: virial definition: 0: Bryan-Norman; 1: 200c; 2: 200b.
% cn: Rv/R_{-2} for subhalos

G=43.0071;
HUBBLE0=100;  %km/s/(Mpc/h)
Omega0=0.3;OmegaLambda=0.7;scaleF=1./(1+z);
Hz=HUBBLE0 * sqrt(Omega0 ./scaleF.^3+ (1 -Omega0 -OmegaLambda) ./ scaleF.^2 +OmegaLambda);

Hratio=Hz/HUBBLE0;
OmegaZ=Omega0./scaleF.^3./Hratio.^2;
switch virtype   %virial factor with respect to critical density
    case 0
        virialF=18.0*pi^2+82.0*(OmegaZ-1)-39.0*(OmegaZ-1).^2;
        A=7.85;
	B=-0.081;
	C=-0.71;
    case 1
        virialF=200.;
        A=5.71;
	B=-0.084;
	C=-0.47;
    case 2'
        A=10.14;
	B=-0.081;
	C=-1.01;
        virialF=200*OmegaZ;
    otherwise
        error('virialtype must be 0/1/2')
end

Mp=2e2; %10^10Msun/h
c=A*(M/Mp).^B.*(1+z).^C;

rhoc=(3*Hz.^2)/(8*pi*G);
rhos=virialF/3.*c.^3./(log(1+c)-c./(1+c))/OmegaZ; %in units of rhob
rv=(M./(4*pi/3*virialF.*rhoc)).^(1/3);

sig=rhos./c./(1+c).^2.*exp(2*cn)*2.*r.*besselk(1,2*r./(rv./cn));