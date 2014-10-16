ns=1;
OmegaM=0.3;
h=0.7;
fb=0.168;
sig8=0.8;
z=0;

wx=@(x) 3*(sin(x)-x.*cos(x))./x.^3;
delt8=quadgk(@(k) wx(k*8).^2.*k.^(2+ns).*TF_BBKS(OmegaM,h,fb,k).^2,0,inf,'RelTol',1e-2);  %this integral is contributed primarily by k from 0.1 to 1
Ap=sig8.^2/delt8*pi.^2*2; %normalization

fun=@(k) Ap*TF_BBKS(OmegaM,h,fb,k).^2.*k.^ns;
% fun=@(x) exp(-x.^2/2)*pi*2;
% fun2=@(x) exp(-x.^2*pi);

T=1;
N=2048;
dk=T/N;
[kx,ky]=meshgrid((0:N)*dk,(0:N)*dk);
k=sqrt(kx.^2+ky.^2);
Pk=fun(k);
Pk=[Pk,Pk(:,end-1:-1:2)];
Pk=[Pk;Pk(end-1:-1:2,:)];
Pk(1)=0;

dr=2*pi./(2*N)/dk;
r=dr*(0:N);
[rx,ry]=meshgrid(r,r);
r=sqrt(rx.^2+ry.^2);
ksi=ifft2(Pk/dr/dr,'symmetric');  %2d FT: (1/dr)^2 difference
%%
ksi2=proj_corr_lin(r(1,:));

%%
figure;loglog(r(1,:),abs(ksi(1,1:N+1)),'gx')
hold on;
% loglog(r(1,:),fun2(r(1,:)),'ro')
plot(r(1,:),abs(ksi2),'r.');
loglog(r(1,:),(ksi(1,1:N+1)),'g-')
hold on;
% loglog(r(1,:),fun2(r(1,:)),'ro')
plot(r(1,:),(ksi2),'r-');

xlabel('r_p/(Mpc/h)');ylabel('w');
print('-depsc','proj_corr_lin_fft.eps');

%%
figure;loglog(r(1,:),(ksi(1,1:N+1)),'.')
hold on;
% loglog(r(1,:),fun2(r(1,:)),'ro')
plot(r(1,:),(ksi2),'ro');

