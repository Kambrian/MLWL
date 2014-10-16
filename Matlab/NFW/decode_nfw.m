function [rhos,rs,rv]=decode_nfw(M,c,z,virtype)
% DeltaSigma=sigmacrit*shear for NFW halo,physical value, ref: Wright & Brainerd 2000
%
% M: 10^10Msun/h
% c: concentration
% z: redshift
% virtype: virial definition: 0: Bryan-Norman; 1: 200c; 2: 200b.
% return rhos,rs,rv
% rv: physical virial radius, Mpc/h

G=43.0071;
HUBBLE0=100;  %km/s/(Mpc/h)
Omega0=0.3;OmegaLambda=0.7;scaleF=1./(1+z);
Hz=HUBBLE0 * sqrt(Omega0 ./scaleF.^3+ (1 -Omega0 -OmegaLambda) ./ scaleF.^2 +OmegaLambda);

Hratio=Hz/HUBBLE0;
OmegaZ=Omega0./scaleF.^3./Hratio.^2;
switch virtype   %virial factor with respect to critical density
    case 0
        virialF=18.0*pi^2+82.0*(OmegaZ-1)-39.0*(OmegaZ-1).^2;
    case 1
        virialF=200.;
    case 2
        virialF=200*OmegaZ;
    otherwise
        error('virialtype must be 0/1/2')
end
	
rhoc=(3*Hz.^2)/(8*pi*G);
rhos=virialF/3.*c.^3./(log(1+c)-c./(1+c)).*rhoc;
rv=(M./(4*pi/3*virialF.*rhoc)).^(1/3);

rs=rv./c;

