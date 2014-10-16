function n=smfunc(logm,flag_new)
% stellar mass function in the local universe
% Li&White 09, stepwise Schechter func fit.
% input: logm=log10(m), where m is stellar mass in Msun/h^2
% output: dN/dV/dlog10M: 1/(Mpc/h)^3

if nargin<2
    flag_new=0;
end
if(~flag_new)
%Li&White09
fai=[0.0146,0.0132,0.0044,0.0083];
alpha=1+[-1.13,-0.90,-1.99,-1.155];
logm0=[9.61,10.37,10.71,10.525];

n=logm;
f=logm<9.33;
n(f)=log(10)*fai(1)*10.^((logm(f)-logm0(1))*alpha(1)).*exp(-10.^(logm(f)-logm0(1)));
f=logm>=9.33&logm<10.67;
n(f)=log(10)*fai(2)*10.^((logm(f)-logm0(2))*alpha(2)).*exp(-10.^(logm(f)-logm0(2)));
f=logm>=10.67;
n(f)=log(10)*fai(3)*10.^((logm(f)-logm0(3))*alpha(3)).*exp(-10.^(logm(f)-logm0(3)));

%overall fit
% n=log(10)*fai(4)*10.^((logm-logm0(4))*alpha(4)).*exp(-10.^(logm-logm0(4)));
else
%Guo et al.,10
logh2=log10(0.73^2);%convert stellar mass unit
fai=[0.0159,0.0121,0.0032];
alpha=1+[-1.11,-0.938,-2.33];
logm0=[9.84,10.71,11.09];

n=logm;
logm=logm-logh2; %from Msun/h^2 to Msun
f=logm<9.60;
n(f)=log(10)*fai(1)*10.^((logm(f)-logm0(1))*alpha(1)).*exp(-10.^(logm(f)-logm0(1)));
f=logm>=9.60&logm<10.94;
n(f)=log(10)*fai(2)*10.^((logm(f)-logm0(2))*alpha(2)).*exp(-10.^(logm(f)-logm0(2)));
f=logm>=10.94;
n(f)=log(10)*fai(3)*10.^((logm(f)-logm0(3))*alpha(3)).*exp(-10.^(logm(f)-logm0(3)));
% n=log10(0.0032*(1e12/10^11.09).^(1-2.33).*exp(-1e12/10.^11.09));
end