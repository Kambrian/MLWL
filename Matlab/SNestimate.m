zl=0.2;zs=0.2:0.5:1000;
figure;
plot(zs,sigma_crit(0.3,0.7,zl,zs).^-2);

f=@(z) sigma_crit(0.3,0.7,0.2,z).^-2.*source_zdistr(z);
nf=@(z) f(z)/f(0.3);
g=@(z) source_zdistr(z);
quad(nf,0.2,10)/quad(g,0.2,10)*f(0.3)

%%
z=0.1:0.1:0.5;s=z;ss=s;
for i=1:numel(z)
    s(i)=sigma_crit_eff(z(i));
    ss(i)=sigma_crit_eff_normal(z(i));
end
myfigure;plot(z,s,'r',z,ss,'b');
xlabel('$z_{s}$');
ylabel('$\Sigma_{c,eff}$');
legend('sigma weighted','usual');
print('-depsc','sigma_crit_eff.eps');

d=comoving_function(comoving_dist(0.3,0.7,z),1)./(1+z);
myfigure;plot(z,s.*d,'r',z,ss.*d,'b');
xlabel('$z_{s}$');
ylabel('$D_l\cdot\Sigma_{c,eff}$');
legend('sigma weighted','usual');
print('-depsc','sigma_crit_effxDl.eps');

%%
M=1e4;zl=0.3;
n0=1.2; %source number density, arcmin^-2;
nlens=11000/3;
nstack=nlens*n0; %stacked number density
errchi=0.4*sqrt(2); % noise of chi ellipticity, including intrinsic and measurement, two component value
rmin=0.2;rmax=2;

[sigav,r200]=nfw_avDeltSig([rmin,rmax],M,zl);
dl=comoving_function(comoving_dist(0.3,0.7,zl),1)/(1+zl);
r=logspace(log10(rmin*r200),log10(rmax*r200),20);  % logscale radial bin
% has good over S/N trend
% r=linspace((rmin*r200),(rmax*r200),20); % linear scale radial has highest inner S/N
sig=nfw_DeltSig(r,M,zl,0);
theta=r/dl/pi*180*60; %angular seperation in arcminutes
sn=sig(1:end-1).*sqrt(pi*diff(theta.^2))*sqrt(nstack)*2*sqrt(2)/sigma_crit_eff(zl)/errchi;

[snav,snav_n,snmass]=SNmass(zl,M,nlens);
myfigure;semilogx(r(2:end),sn);hold on; 
plot(r,repmat(snav,size(r)),'r-');
plot(r,repmat(snav_n,size(r)),'g-');
plot(r,repmat(snmass,size(r)),'k-');
l=legend('$\Delta\Sigma_{wght}(r)$','$<\Delta\Sigma_{wght}>$','$<\Delta\Sigma>$','$M_{\zeta}$');
set(l,'interpreter','latex','location','west','box','off');
xlabel('R/(Mpc/h)');ylabel('S/N');
%print('-depsc','SN_single.eps');
%%
%%
nbin=[10,3];
f=cat.Mass>0;
f=f&cat.Mult>=2;
[xx,yy,n,s]=densitygrid(log10(cat.Mass(f)),cat.MedianZ(f),nbin);
myfigure;
C=contour(xx,yy,n,[1,10,50,100,300,500,800,1000,1500,2000,3000,4000,5000,8000]);
clabel(C,'fontsize',15);
xlabel('$log(M_{dyn}/(M_{\odot}/h))$');
ylabel('$z$');
% title(['Ngroups, ',num2str(nbin),'$\times$',num2str(nbin),' grid']);
print('-depsc',['nbin',num2str(nbin(1)),'x',num2str(nbin(2)),'.eps']);

sn=n;
sn_n=n;snmass=n;
for i=1:nbin(2)
    for j=1:nbin(1)
     [sn(i,j),sn_n(i,j),snmass(i,j)]=SNmass(yy(i,j),10^(xx(i,j)-10),n(i,j));
    end
end
myfigure;
[C,h]=contour(xx,yy,sn,[1,5,10,15,20,30]);
clabel(C,'fontsize',15);
% imagesc(minmax(catmock.VelDisp'),minmax(catmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
% colorbar;
l=legend('$S/N,<\Delta\Sigma_{wght}>$');set(l,'interpreter','latex','location','northwest');
xlabel('$log(M_{dyn}/(M_{\odot}/h))$');
ylabel('$z$');
% print('-depsc',['sn_cont',num2str(nbin(1)),'x',num2str(2),'.eps']);

myfigure;
[C,h]=contour(xx,yy,sn_n,[1,5,10,15,20,30]);
clabel(C,'fontsize',15);
% imagesc(minmax(catmock.VelDisp'),minmax(catmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
% colorbar;
l=legend('$S/N,<\Delta\Sigma>$');set(l,'interpreter','latex','location','northwest');
xlabel('$log(M_{dyn}/(M_{\odot}/h))$');
ylabel('$z$');
% print('-depsc',['sn_cont_norm',num2str(nbin(1)),'x',num2str(2),'.eps']);

myfigure;
[C,h]=contour(xx,yy,snmass,[1,3,5,8,12,15,20,30]);
clabel(C,'fontsize',15);
% imagesc(minmax(catmock.VelDisp'),minmax(catmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
% colorbar;
l=legend('$S/N,M_{\zeta}$');set(l,'interpreter','latex','location','northwest');
xlabel('$log(M_{dyn}/(M_{\odot}/h))$');
ylabel('$z$');
print('-depsc',['sn_cont_mass',num2str(nbin(1)),'x',num2str(nbin(2)),'.eps']);

%%
% cat.Mass=sqrt(cat.Mass.*cat.LumMass);
mbin=7;
nbin=[mbin,3];
f=cat.Mass>0;
f=f&cat.Mult>=2;
[xx,yy,n,s]=densitygrid(log10(cat.Mass(f)),cat.MedianZ(f),nbin);
myfigure;
C=contour(xx,yy,n,[1,10,50,100,300,500,800,1000,1500,2000,3000,4000,5000,8000]);
clabel(C,'fontsize',15);
xlabel('$log(M_{dyn}/(M_{\odot}/h))$');
ylabel('$z$');
% title(['Ngroups, ',num2str(nbin),'$\times$',num2str(nbin),' grid']);
% print('-depsc',['nbin',num2str(nbin(1)),'x',num2str(nbin(2)),'.eps']);

sn=n;
sn_n=n;snmass=n;
for i=1:nbin(2)
    for j=1:nbin(1)
     [sn(i,j),sn_n(i,j),snmass(i,j)]=SNmass(yy(i,j),10^(xx(i,j)-10),n(i,j));
    end
end

myfigure;
semilogy(xx(3,:),snmass(3,:),'ro-','displayname',['z=',num2str(yy(3,1),'%2.1f')]);
hold on;
semilogy(xx(2,:),snmass(2,:),'go-','displayname',['z=',num2str(yy(2,1),'%2.1f')]);
semilogy(xx(1,:),snmass(1,:),'bo-','displayname',['z=',num2str(yy(1,1),'%2.1f')]);

zl=0.2;
nbin=mbin+1;
f=cat.Mass>0;
f=f&cat.Mult>=2;
[x,n]=loghist(cat.Mass(f),nbin);
sn=n;sn_n=n;snmass=n;
for i=1:nbin-1
    [sn(i),sn_n(i),snmass(i)]=SNmass(zl,x(i)/1e10,n(i));
end
% myfigure;
plot(log10(x),snmass,'ko-','displayname','all at z=0.2');
l=legend('show');
set(l,'interpreter','latex','location','southeast','box','off');
xlabel('$log(M_{dyn}/(M_{\odot}/h))$');
ylabel('$SN_{M_{\zeta}}$');
% print('-depsc',['sn_mass',num2str(mbin),'z.eps']);
%%
nbin=[10,10];
f=cat.Mass>0;
f=f&cat.Mult>=2;
[xx,yy,n,s]=densitygrid(log10(cat.Mass(f)),cat.MedianZ(f),nbin);
myfigure;
C=contour(xx,yy,n,[1,10,50,100,300,500,800]);
clabel(C,'fontsize',15);
xlabel('$log(M_{dyn}/(M_{\odot}/h))$');
ylabel('$z$');
% print('-depsc',['nbin',num2str(nbin(1)),'x',num2str(nbin(2)),'.eps']);

sn=n;
sn_n=n;snmass=n;
for i=1:nbin(2)
    for j=1:nbin(1)
     [sn(i,j),sn_n(i,j),snmass(i,j)]=SNmass(yy(i,j),10^(xx(i,j)-10),sum(n(:)));
    end
end
myfigure;
[C,h]=contour(xx,yy,snmass,[1,3,10,30,100,300]);
clabel(C,'fontsize',15);
l=legend('$S/N,M_{\zeta}$');set(l,'interpreter','latex','location','northwest');
xlabel('$log(M_{dyn}/(M_{\odot}/h))$');
ylabel('$z$');
print('-depsc',['sn_cont_mass_Mad',num2str(nbin(1)),'x',num2str(nbin(2)),'.eps']);