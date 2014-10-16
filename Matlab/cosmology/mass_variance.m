function sigm2=mass_variance(m,OmegaM,h,fb,sig8)
% sigma(m)^2, mass variance in linear power spectrum
% m: in units 10^10M_sun/h


G=43.0071;
HUBBLE0=100;  %km/s/(Mpc/h)
rhom=OmegaM*3*HUBBLE0^2/(8*pi*G);
Rh=(3*m/4/pi/rhom).^(1/3); %filter radius
ns=1;  %primordial index

wx=@(x) 3*(sin(x)-x.*cos(x))./x.^3;
sigm2=m;
for i=1:numel(m)
sigm2(i)=quadgk(@(k) wx(k*Rh(i)).^2.*k.^(2+ns).*TF_BBKS(OmegaM,h,fb,k).^2,0,inf,'RelTol',1e-2);
end

delt8=quadgk(@(k) wx(k*8).^2.*k.^(2+ns).*TF_BBKS(OmegaM,h,fb,k).^2,0,inf,'RelTol',1e-2);  %this integral is contributed primarily by k from 0.1 to 1
sigm2=sigm2/delt8*sig8.^2;