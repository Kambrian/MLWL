function [sig,r200]=nfw_avDeltSig(r,M,z)
% average <DeltaSigma>=<sigmacrit*shear> for NFW halo,physical value, ref: Wright & Brainerd 2000
% r: [rmin,rmax], in units of rvir
% M: 10^10Msun/h
% z: redshift

G=43.0071;
HUBBLE0=100;  %km/s/(Mpc/h)
Omega0=0.3;OmegaLambda=0.7;scaleF=1./(1+z);
Hz=HUBBLE0 * sqrt(Omega0 ./scaleF.^3+ (1 -Omega0 -OmegaLambda) ./ scaleF.^2 +OmegaLambda);


Mc=1.5e3; %M_{*},in 10^10Msun/h, characteristic mass scale
c=9./(1+z).*(M/Mc).^-0.13;
rhoc=(3*Hz.^2)/(8*pi*G);
rhos=200/3*c.^3./(log(1+c)-c./(1+c)).*rhoc;
r200=(M./(4*pi/3*200*rhoc)).^(1/3);
rs=r200./c;

rmin=r(1)*c;
rmax=r(2)*c;
gl=@(x) 8*atanh(sqrt((1-x)./(1+x)))./x.^2./sqrt(1-x.^2)+4./x.^2.*log(x/2)...
        -2./(x.^2-1)+4*atanh(sqrt((1-x)./(1+x)))./(x.^2-1)./sqrt(1-x.^2);
gr=@(x) 8*atan(sqrt((x-1)./(1+x)))./x.^2./sqrt(x.^2-1)+4./x.^2.*log(x/2)...
        -2./(x.^2-1)+4*atan(sqrt((x-1)./(1+x)))./(x.^2-1).^(3/2);   
if rmin<1
    if rmax>1
    y=quadgk(@(x) gl(x).*2*pi.*x,rmin,1)+quadgk(@(x) gr(x)*2*pi.*x,1,rmax);
    else
    y=quad(@(x) gl(x)*2*pi.*x,rmin,rmax);
    end
else
    y=quad(@(x) gr(x).*2*pi.*x,rmin,rmax);
end
    
sig=y/pi/(rmax^2-rmin^2)*rs*rhos;


% loglog(rx,sig);
