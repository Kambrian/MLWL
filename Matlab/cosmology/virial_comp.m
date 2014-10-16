function [M,R,virF]=virial_comp(z,Omega0)
% scaling relation between different virial definitions when applied on the
% same halo
% M~rho*r^3, for [vir,200c,200b]
% returned M and R has arbitrary normalization, but the correct scaling
% between different definitions.
%
% Note the scaling with redshift for the same definition is not
% guanranteed.
%
% Ref:arxiv:0911.0436 appendix. Giocoli et al.

if nargin<2
    Omega0=0.3;
end
OmegaLambda=1-Omega0;
scaleF=1./(1+z);
% G=43007.1;
HUBBLE0=0.1;

Hz=HUBBLE0 * sqrt(Omega0 ./scaleF.^3+ (1 -Omega0 -OmegaLambda) ./ scaleF.^2 +OmegaLambda);
Hratio=Hz/HUBBLE0;
OmegaZ=Omega0./scaleF.^3./Hratio.^2;
virialF=18.0*pi^2+82.0*(OmegaZ-1)-39.0*(OmegaZ-1).^2;
virialF0=18.0*pi^2+82.0*(Omega0-1)-39.0*(Omega0-1).^2;

Rv=ones(size(z));
R200c=0.746*(virialF./virialF0).^0.395.*Rv;
R200b=1.236*(virialF./virialF0).^-0.438.*Rv;

Mv=virialF.*Rv.^3;
M200c=200*R200c.^3;
M200b=200*OmegaZ.*R200b.^3;

M=[Mv,M200c,M200b];
R=[Rv,R200c,R200b];
virF=[virialF,200,200*OmegaZ];