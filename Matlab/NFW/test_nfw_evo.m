dir='/home/kam/Projects/Lensing/data';
z=0:0.1:1;
M=1e3;

[rhos0,rs0,rv0]=nfw_evo(M,z,0);
[rhos1,rs1,rv1]=nfw_evo(M,z,1);
[rhos2,rs2,rv2]=nfw_evo(M,z,2);
c0=rv0./rs0;c1=rv1./rs1;c2=rv2./rs2;
M0=4*pi*rs0.^3.*rhos0.*(log(1+c0)-c0./(1+c0));
M1=4*pi*rs1.^3.*rhos1.*(log(1+c1)-c1./(1+c1));
M2=4*pi*rs2.^3.*rhos2.*(log(1+c2)-c2./(1+c2));
myfigure;
plot(z,rhos0,'r-o');
hold on;
plot(z,rhos1,'g-o');
plot(z,rhos2,'b-o');
plot(z,rhos0./(1+z).^3,'r--o');
plot(z,rhos1./(1+z).^3,'g--o');
plot(z,rhos2./(1+z).^3,'b--o');
xlabel('z');ylabel('$\rho_s$');
l=legend('$\Delta_{th}$','200c','200b');
set(l,'interpreter','latex');
print('-depsc',[dir,'/nfw_rhos_evo.eps']);

myfigure;
plot(z,rs0,'r-o');
hold on;
plot(z,rs1,'g-o');
plot(z,rs2,'b-o');
plot(z,rs0.*(1+z),'r--o');
plot(z,rs1.*(1+z),'g--o');
plot(z,rs2.*(1+z),'b--o');
xlabel('z');ylabel('$r_s$');
l=legend('$\Delta_{th}$','200c','200b');
set(l,'interpreter','latex');
print('-depsc',[dir,'/nfw_rs_evo.eps']);

myfigure;
plot(z,rv0,'r-o');
hold on;
plot(z,rv1,'g-o');
plot(z,rv2,'b-o');
plot(z,rv0.*(1+z),'r--o');
plot(z,rv1.*(1+z),'g--o');
plot(z,rv2.*(1+z),'b--o');
xlabel('z');ylabel('$r_v$');
l=legend('$\Delta_{th}$','200c','200b');
set(l,'interpreter','latex');
print('-depsc',[dir,'/nfw_rv_evo.eps']);

f=@(x) log(1+x)-x./(1+x);
x=c2;
figure;plot(z,f(x),'o-');
figure;plot(z,rs2.^3.*f(x),'o-');
%%
nfw=@(rhos,x) rhos./x./(1+x).^2;
r=logspace(-4,1,10);
myfigure;
colors=colormap(jet(numel(z)));
% physical density in physical radius
for i=1:3:numel(z)
    loglog(r,r.^1.*nfw(rhos0(i),r./rs0(i)),'--','color',colors(i,:));
    hold on;
end
% comoving density in comoving radius
for i=1:3:numel(z)
%     loglog(r*(1+z(i)),nfw(rhos0(i)/(1+z(i))^3,r./rs0(i)),'color',colors(i,:));
    loglog(r,r.^1.*nfw(rhos0(i)/(1+z(i))^3,r./(rs0(i)*(1+z(i)))),'o','color',colors(i,:));
    hold on;
end
plot([rs0(1),rs0(1)],[1,10]);
xlabel('r/(Mpc/h)');
ylabel('$r\rho$');
print('-depsc',[dir,'/nfw_rho_evo.eps']);
%%
x=logspace(-3,1,10);
for i=1:3:numel(z)
    r=x*rv0(i);
    loglog(x,r.^1.*nfw(rhos0(i),r./rs0(i)),'--','color',colors(i,:));
    hold on;
end
% comoving density in comoving radius
for i=1:3:numel(z)
    r=x*rv0(i)*(1+z(i)); %comoving coord
%     loglog(r*(1+z(i)),nfw(rhos0(i)/(1+z(i))^3,r./rs0(i)),'color',colors(i,:));
    loglog(x,r.^1.*nfw(rhos0(i)/(1+z(i))^3,r./(rs0(i)*(1+z(i)))),'o','color',colors(i,:));
    hold on;
end

xlabel('r/rv');
ylabel('$r\rho$');

%%
M=1e3;
r=logspace(-4,1,20);
myfigure;
colors=colormap(jet(numel(z)));
% physical density in physical radius
for i=1:3:numel(z)
    loglog(r,nfw_DeltSig(r,M,z(i),0),'--','color',colors(i,:));
    hold on;
end
% comoving density in comoving radius
for i=1:3:numel(z)
    loglog(r,nfw_DeltSig(r./(1+z(i)),M,z(i),0)/(1+z(i))^2,'o','color',colors(i,:));
    hold on;
end
% plot([rs0(1),rs0(1)],[1,10]);
xlabel('r/(Mpc/h)');
ylabel('$\Delta\Sigma$');
print('-depsc',[dir,'/nfw_surfdens_evo.eps']);

%% similar sensitivity to halo mass
M=logspace(2,4,10);
z=0.1;
r=logspace(-4,1,20);
myfigure;
colors=colormap(jet(numel(M)));
% physical density in physical radius
for i=1:numel(M)
    loglog(r,nfw_DeltSig(r,M(i),z,0),'--','color',colors(i,:),'displayname',['$',printexp10(M(i)*1e10),'M\odot/h','$']);
    hold on;
end
l=legend('show');set(l,'interpreter','latex','location','southwest');
% comoving density in comoving radius
for i=1:numel(M)
    loglog(r,nfw_DeltSig(r./(1+z),M(i),z,0)/(1+z)^2,'o','color',colors(i,:));
    hold on;
end

% plot([rs0(1),rs0(1)],[1,10]);
xlabel('r/(Mpc/h)');
ylabel('$\Delta\Sigma$');
print('-depsc',[dir,'/nfw_surfdens_sensitivity.eps']);
%%
r=[100,1000,2000]/1000;
colors='rgbkc';
M=logspace(2,4,100);
% r=comoving_200b_radius(M,0.3);
z=0.2;
dsig=[];
myfigure;
for i=1:numel(r)
    s=nfw_DeltSig(r(i),M,z,2);
    dsig=[dsig;s];
    ds=diff(log(M))./diff(log(s));
    semilogx(M(1:end-1)*1e10,ds,'-','color',colors(i));
%     plot(log(M),log(s),'-','color',colors(i));
    hold on;
end
l=legend('0.1','1','2'); set(l,'location','northwest','interpreter','latex');
% ylim([0.98,1.1]);
xlabel('$M[\rm{M}_\odot/h]$');
ylabel('$\rm{d}\ln M/\rm{d}\ln\Delta\Sigma$');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/MassResponse.eps');
%%
colors='rgbkc';
M=logspace(2,4,100);
r=comoving_200b_radius(M,0.3);
z=0.2;
dsig=[];
myfigure;
for i=1:numel(r)
    s=nfw_surf_overdensity(r(i),M(i),z,2);
    dsig=[dsig;s];
%     ds=diff(log(M))./diff(log(s));
%     semilogx(M(1:end-1)*1e10,ds,'-','color',colors(i));
%     plot(log(M),log(s),'-','color',colors(i));
%     hold on;
end
ds=diff(log(M))./diff(log(dsig'));
plot(log10(M(1:end-1)),ds,'-');




