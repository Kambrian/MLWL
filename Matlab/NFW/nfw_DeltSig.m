function [sig,rv]=nfw_DeltSig(r,M,z,virtype,c)
% DeltaSigma=sigmacrit*shear for NFW halo,physical value, ref: Wright & Brainerd 2000
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
    case 2
        virialF=200*OmegaZ;
        A=10.14;
	B=-0.081;
	C=-1.01;
    otherwise
        error('virialtype must be 0/1/2')
end
	
% Mc=1.5e3; %M_{*},in 10^10Msun/h, characteristic mass scale at z=0
%Mc=charact_mass(z,0.3,0.7,0.16,0.8); %M_{*},in 10^10Msun/h, characteristic mass scale at z
%c=9./(1+z).*(M/Mc).^-0.13;  %% this should be different for diff virtypes=========================== to be updated

if nargin<5||c<=0 %only set c if it's not specified or c<=0
% Duffy et al., 2008
Mc=2e2; %M_{*},in 10^10Msun/h, characteristic mass scale at z
c=A.*(1+z).^C.*(M/Mc).^B;  %% this should be different for diff virtypes=========================== to be updated
%c=c/3;
end

rhoc=(3*Hz.^2)/(8*pi*G);
rhos=virialF/3.*c.^3./(log(1+c)-c./(1+c)).*rhoc;
rv=(M./(4*pi/3*virialF.*rhoc)).^(1/3);

rs=rv./c;

rx=r./rs;
sig=zeros(size(rx));

gl=@(x) 8*atanh(sqrt((1-x)./(1+x)))./x.^2./sqrt(1-x.^2)+4./x.^2.*log(x/2)...
        -2./(x.^2-1)+4*atanh(sqrt((1-x)./(1+x)))./(x.^2-1)./sqrt(1-x.^2);
gr=@(x) 8*atan(sqrt((x-1)./(1+x)))./x.^2./sqrt(x.^2-1)+4./x.^2.*log(x/2)...
        -2./(x.^2-1)+4*atan(sqrt((x-1)./(1+x)))./(x.^2-1).^(3/2);   

    sig(rx<1)=gl(rx(rx<1));
    sig(rx==1)=10/3+4*log(1/2);
%     sig(rx>1&rx<5*c)=gr(rx(rx>1&rx<5*c));
    sig(rx>1)=gr(rx(rx>1));
    
%     for i=1:numel(rx)
%         x=rx(i);
%         if x<1
%             sig(i)=gl(x);
%         else if x==1
%                 sig(i)=10/3+4*log(1/2);
%             else if x<10*c
%                 sig(i)=gr(x);
%                 else
%                     sig(i)=0;
%                 end
%             end
%         end
%     end
sig=sig.*rs.*rhos;


% loglog(rx,sig);
