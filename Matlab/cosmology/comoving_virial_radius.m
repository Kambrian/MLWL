function Rvir=comoving_virial_radius(Mvir,scaleF,vir_type)
%function Rvir=comoving_virial_radius(Mvir,scaleF,vir_type)
%Note input Mvir is physical rather than particle number, in units
%10^10Msun/h
% output: Rvir, unit: kpc/h
%Omega0,OmegaLambda being density param at z=0
%         scaleF=a(snapnum+1);
%         vir_type: optional
%                       0: default, tophat (Bryan&Norman)
%                       1: 200c
%                       2: 200b
global Omega0 OmegaLambda
if isempty(Omega0)
    Omega0=0.3;
    disp('Cosmology: Omega0=0.3');
end
if isempty(OmegaLambda)
    OmegaLambda=1-Omega0;
end
G=43007.1;
HUBBLE0=0.1;
Hz=HUBBLE0 * sqrt(Omega0 ./scaleF.^3+ (1 -Omega0 -OmegaLambda) ./ scaleF.^2 +OmegaLambda);
Hratio=Hz/HUBBLE0;
OmegaZ=Omega0./scaleF.^3./Hratio.^2;

if nargin<3
    vir_type=0;
end

switch vir_type
    case 0
        virialF=18.0*pi^2+82.0*(OmegaZ-1)-39.0*(OmegaZ-1).^2
    case 1
        virialF=200;
    case 2
        virialF=200*OmegaZ;
    case 3
        virialF=6*OmegaZ;
    otherwise
        error('wrong vir_type');
end

Rvir=(2.0*G*Mvir./virialF./Hz.^2).^(1.0/3)./scaleF;