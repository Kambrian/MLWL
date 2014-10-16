function [sig,rv]=nfw_surf_overdensity(r,M,z,virtype)
%Sigma=int(rho_NFW,z,-inf,inf)/rho_background for NFW halo,
% surface mass overdensity for a NFW halo sunperimposed onto a uniform background
% ref: Wright & Brainerd 2000
% rv: physical virial radius
%
% r: Mpc/h, physical
% M: 10^10Msun/h
% z: redshift
% virtype: virial definition: 0: Bryan-Norman; 1: 200c; 2: 200b.

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

% Mc=1.5e3; %M_{*},in 10^10Msun/h, characteristic mass scale at z=0
%Mc=charact_mass(z,0.3,0.7,0.154,0.7); %M_{*},in 10^10Msun/h, characteristic mass scale at z
%c=9./(1+z).*(M/Mc).^-0.13;  %% this should be different for diff virtypes=========================== to be updated
rhoc=(3*Hz.^2)/(8*pi*G);
rhos=virialF/3.*c.^3./(log(1+c)-c./(1+c))/OmegaZ;
rv=(M./(4*pi/3*virialF.*rhoc)).^(1/3);

rs=rv./c;

rx=r./rs;
sig=zeros(size(rx))+2./3;

fl=@(x) 2./(x.^2-1).*(1-2./sqrt(1-x.^2).*atanh(sqrt((1-x)./(1+x))));
fr=@(x) 2./(x.^2-1).*(1-2./sqrt(x.^2-1).*atan(sqrt((x-1)./(1+x))));

sig(rx<1)=fl(rx(rx<1));
sig(rx>1)=fr(rx(rx>1));
%    for i=1:numel(rx)
%         x=rx(i);
%         if x<1
%             sig(i)=fl(x);
%         else if x==1
%                 sig(i)=2./3;
%             else
%                 sig(i)=fr(x);
%             end
%         end
%     end
sig=sig.*rs.*rhos;


% loglog(rx,sig);
