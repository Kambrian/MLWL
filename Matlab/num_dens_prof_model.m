H=0.1;
c=3e5;
sig=0.1*c/H;
f=@(r,z,zp) 1./(r.^2+z.^2).*exp(-(zp-z).^2/2/sig.^2);
r=logspace(0,5,20);
y=r;
for i=1:20
    y(i)=quad2d(@(z,zp) f(r(i),z,zp),0,1e8,0,1e8);
end
%%
macro=default_WLparam;
macro.ZGAMA=0;
GAMA_WL_init;
z=[e{1}(:,6);e{2}(:,6);e{3}(:,6)];
zbrd=[e{1}(:,7:8);e{2}(:,7:8);e{3}(:,7:8)];
catall=structcat([grp{1},grp{2},grp{3}]);
f=catall.Mult>=macros.MULT_MIN;
zmean=zeros(4,1);halomass_mean=zmean;halomass_err=zmean;
for i=1:4
zmean(i)=mean(catall.Zfof(f&catall.mid==i));
halomass_mean(i)=mean(catall.HaloMass(f&catall.mid==i));
halomass_err(i)=std(catall.HaloMass(f&catall.mid==i));
end
M=halomass_mean/1e10;eM=halomass_err/1e10;
%%
% zerr=abs(diff(zbrd,1,2))/2;
zerr=sqrt(mean((zbrd-repmat(z,1,2)).^2,2));
[xx,yy,n]=densitygrid(z,zerr,[30,30]);
[zmed,zerrmed,~,~,~,~,count]=skeleton(z,zerr,30,0.683);
figure;
contourf(xx,yy,log(n+1));
hold on;
x=0:0.1:1;
plot(x,0.1*x+0.06,'k','linewidth',3);
plot(zmed,zerrmed,'y');
xlabel('z');ylabel('\sigma_z');
% title('CC2');
title('ZEBRA');
legend('ln(n+1)','\sigma_z=0.1*(z+0.6)','median');
% print('-depsc','zerr_z_ZEBRA.eps');
% figure;
% plot(z,zerr,'.');
% figure;hist(zerr./z,0:0.05:1);xlim([0,1]);
% figure;hist(zerr./(1+z),0:0.05:1);xlim([0,1]);
% figure;hist(zerr);

% hold on;
% plot(z,zerr(:,2),'.');
%%
ltyp={'ro-','ro--';'g>-','g>--';'bs-','bs--';'kd-','kd--';};
ldtyp={'ro:';'g>:';'bs:';'kd:';};
lltyp={'ro-';'g>-';'bs-';'kd-';};
mtyp={'ro';'g>';'bs';'kd';};
color={'r';'g';'b';'k'};
%%
zbin=linspace(min(z),max(z),50);
n=histc(z,zbin);
v=AD_dist_flat(0.3,0,zbin).^3*(4/180*pi)*(12/180*pi)/3.*(1+zbin).^3; %comoving volume
dndv=n(1:end-1)'./diff(v);
zm=(zbin(1:end-1)+zbin(2:end))/2;
figure;
plot(zm,log(dndv)-log(dndv(1)),'o');
cftool(zm,log(dndv))
k=mean(log(dndv/dndv(1))./zm)
% k=sum(log(dndv/dndv(1)))/sum(zbin(1:end-1))
hold on;
plot(zm,k*zm,'r');
figure;
plot(comoving_dist(0.3,0.7,zbin(1:end-1)),log(dndv));

figure;
plot(zbin,comoving_dist(0.3,0.7,zbin))
hold on;
plot(zbin,3000*zbin)

%%
r=logspace(-1,1,10);
% z=linspace(0,0.5,21);
z=0.1:0.01:0.2;
[rr,zz]=meshgrid(r,z);
y=num_dens_prof_toy(rr,zz,0.16,1e5);
figure;
mesh(log10(rr),zz,y);
%%
r=logspace(-3,1,10);
M=1e4;zc=0.15;
c=3e5;
G=43.0071; % gravity const
H=100; %hubble param, km/s/(Mpc/h)
virialF_b=200; %virial factor, 200 times background
OmegaM=0.3;
rv=(M*2*G/virialF_b/OmegaM/H^2)^(1/3);
req=sqrt(virialF_b/3)*rv;
    
figure;
yp=num_dens_prof_toy2d(r,zc,M,0.1,0.1);
loglog(r,yp-1,'o');

f=@(x) (2./x.*atan(sqrt(1-x.^2)./x)-2*sqrt(1-x.^2))*req;
% f=@(x) pi*req./x;
hold on;
h=@(z) H*sqrt(OmegaM .*(1+z).^3+(1-OmegaM));
plot(r,h(zc+0.1)/c*1.0*f(r/req),'r-');
plot(r,h(zc+0.1)/c*1.0*nfw_surf_overdensity(r/(1+zc),M,zc,2),'k--');
legend('numerical model','analytical approximation','NFW');
title('M=1e14Msun/h, zc=0.15');
xlabel('r/Mpc');ylabel('n/n_{rand}-1')
% print('-depsc','num_dens_prof_model_anal_vs_num.eps');
%%
sigmz=0.1;
zc=0.1;
H=100;OmegaM=0.3;c=3e5;
h=@(z) H*sqrt(OmegaM .*(1+z).^3+(1-OmegaM));
g=@(z,zp) exp(-(z-zp).^2./2/sigmz^2)/sqrt(2*pi)/sigmz;
d=@(z) exp(-8.8*z);
p=@(z) g(z,zc).*d(z);
quad(p,zc+0.1,1.5)/quad(d,zc+0.1,1.5)%/quad(@(z) g(z,0),0,1.5)

OmegaM=0.3;H=100;c=3e5;
z=0:0.01:1.5;
Hz=H * sqrt(OmegaM .*(1+z).^3+(1-OmegaM));
figure;plot(z,c./Hz)
hold on;
plot(z,g(z,0)*3000/g(0,0),'r--');
sigmh=5.98;
zh=-19.9;
A=1.16e7;
RH=A/sqrt(2*pi)/sigmh*exp(-(z-zh).^2/2/sigmh^2);
plot(z,RH,'.');
figure;plot(z,RH.*Hz/c-1)

figure;
semilogy(z,g(z,0),'.')
%%
figure;
yp=cell(4,1);
for i=1:4
yp{i}=num_dens_prof_toy2d(r,zmean(i),M(i),0.1);
semilogx(r,yp{i},color{i});
hold on;
end
%%
r=logspace(-1,1,10);
z=0.1:0.1:0.4;
M=logspace(2,4,3);
myfigure;
yp=cell(4,1);
for i=1:4
    yp{i}=num_dens_prof_toy2d(r,z(i),M(1),0.1);
h1=semilogx(r,yp{i},'s-','color',color{i});
hold on;
end
for i=1:4
    yp{i}=num_dens_prof_toy2d(r,z(i),M(2),0.1);
h2=semilogx(r,yp{i},'o--','color',color{i});
hold on;
end
for i=1:4
    yp{i}=num_dens_prof_toy2d(r,z(i),M(3),0.1);
h3=semilogx(r,yp{i},'x:','color',color{i});
hold on;
end
legend([h1,h2,h3],'M=1e12','1e13','1e14');
xlabel('r/(Mpc/h)');
ylabel('$n/n_b$');
% print('-depsc','nprof_model.eps');
%%
H=100;OmegaM=0.3;c=3e5;virialF_b=200; %virial factor, 200 times background
h=@(z) H*sqrt(OmegaM .*(1+z).^3+(1-OmegaM));
sigmz=0.1;
r=0.1;
z=0.1:0.05:0.5;
M=logspace(3,5,10);
[zz,MM]=meshgrid(z,M);
yp=zz;
ya=yp;%analytical
   
f=@(x) (2./x.*atan(sqrt(1-x.^2)./x)-2*sqrt(1-x.^2));
% f=@(x) pi./x;
for i=1:numel(zz)
    yp(i)=num_dens_prof_toy2d(r,zz(i),MM(i),0.1,sigmz);
    rv=(MM(i)*2*G/virialF_b/OmegaM/H^2)^(1/3);
    req=sqrt(virialF_b/3)*rv;
    ya(i)=1.0*h(zz(i))/c*f(r/req)*req+1;
end
%%
% figure;mesh(zz,log10(MM),yp);
lv=[1.2,1.5,2:6];
myfigure;
[C,h]=contour(zz,log10(MM)+10,yp,lv);
clabel(C,h);
xlabel('z');
ylabel('$log(M/(M_\odot/h))$');
title('$n/n_b(100kpc)$')
% print('-depsc','nprof_model_distr.eps');
myfigure;
[C,h]=contour(zz,log10(MM)+10,ya,lv);
clabel(C,h);
xlabel('z');
ylabel('$log(M/(M_\odot/h))$');
title('$n/n_b(100kpc)$')
print('-depsc','nprof_model_distr_anal.eps');
%%
G=43007.1;OmegaM=0.3;
H=100;c=3e5;
h=@(z) H*sqrt(OmegaM .*(1+z).^3+(1-OmegaM));
Rv=(G*M/OmegaM).^(1/3)/1000;
req=sqrt(200/3)*Rv;
ynfw=cell(4,1);
yiso=cell(4,1);
r=logspace(-2,1,10);
% f=@(x) (2./x.*atan(sqrt(1-x.^2)./x)-2*sqrt(1-x.^2));
f=@(x) pi./x;
for i=1:4
% yp{i}=num_dens_prof_toy2d(r,zmean(i),M(i),0.1);
yiso{i}=h(zmean(i))/c*1.0*f(r/req(i))*req(i);
ynfw{i}=h(zmean(i))/c*1.0*nfw_surf_overdensity(r/(1+zmean(i)),M(i),zmean(i),2);
end
%%
myfigure;
for i=1:4
    eyy=sqrt((datar.eyp{i}./datar.y{i}).^2+(data.ey{i}./data.y{i}).^2).*data.y{i}./datar.y{i};
    errorbar(data.rvar{i},data.y{i}./datar.y{i}-1,eyy,mtyp{i});
    hold on;
    plot(r(r<req(i)),yiso{i}(r<req(i)),'--','color',color{i});
    plot(r(r<req(i)),ynfw{i}(r<req(i)),color{i});
end
ylim([0,5]);
set(gca,'xscale','log','yscale','log');
xlabel('r/req');ylabel('n_{GAMA}/n_{rand}-1','interpreter','tex');
legend('data','isothermal','nfw');
print('-depsc',[dir,'n_prof_ksi_model.eps']);

myfigure;
for i=1:4
    eyy=sqrt((datar.eyp{i}./datar.y{i}).^2+(data.ey{i}./data.y{i}).^2).*data.y{i}./datar.y{i};
    errorbar(data.rvar{i},data.y{i}./datar.y{i},eyy,mtyp{i});
    hold on;
    plot(r,yiso{i}+1,'--','color',color{i});
    plot(r,ynfw{i}+1,color{i});
end
ylim([0,5]);
set(gca,'xscale','log');
xlabel('r/req');ylabel('n_{GAMA}/n_{rand}','interpreter','tex');
legend('data','isothermal','nfw');
print('-depsc',[dir,'n_prof_rat_model.eps']);
%%
