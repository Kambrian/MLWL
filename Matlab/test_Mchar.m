OmegaM=0.3;h=0.7;fb=0.154;
Rh=8;
wx=@(x) 3*(sin(x)-x.*cos(x))./x.^3;
wk=@(k) wx(k*Rh);
% wk=@(k) exp(-(k*Rh).^2/2);
% wk=@(k) 1;

k=logspace(-3,2,100);
% figure;plot(k,wx(k));
figure;semilogx(k,wk(k).^2.*k.^4.*TF_BBKS(OmegaM,h,fb,k).^2);

sig8=0.7;
delt=quadgk(@(k) wk(k).^2.*k.^3.*TF_BBKS(OmegaM,h,fb,k).^2,0,inf,'RelTol',1e-2);  %this integral is contributed primarily by k from 0.1 to 1
A=2*pi^2*sig8/delt;

data=load('CAMB/camb_matterpower.dat');
figure;loglog(data(:,1),data(:,2));
hold on;
plot(data(:,1),A*data(:,1).*TF_BBKS(0.3,0.7,0.154,data(:,1)).^2,'r')

m=logspace(0,5,20);
figure;loglog(m,mass_variance(m,0.3,0.7,0.154,0.8));
hold on;loglog(m,mass_variance(m,0.3,0.7,0.154,0.5),'r');

addpath('/home/kam/Projects/HBT/code/v8.4/anal/Matlab/show');
z=0:0.1:1;
[w,delta,D,D2,D3,w3]=collapse_barrier(1./(1+z),0.3);
figure;plot(z,w,'-');

z=0:0.1:1;
ms=charact_mass(z,0.3,0.7,0.154,0.8)

% fzero(@(m) mass_variance(m,0.3,0.7,0.154,0.8)-1.68,2e3);