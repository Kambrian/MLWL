clear
clc
load G3Cv4up8/mockcat_1.mat
load G3Cv4up8/G3Cv4up8.mat
%%
myfigure;
for vol=1:9
    load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
f=grpmock.Mult>2;
[xm,ym,yl]=skeleton(grpmock.Rad50(f),grpmock.Luminosity(f),logspace(-3,1,30),0.683);
h2=loglog(xm,ym,'g','linewidth',1);hold on;
end
plot(xm,yl,'g--','linewidth',1);
f=grp.Mult>2;
[xm,ym,yl]=skeleton(grp.Rad50(f),grp.Luminosity(f),logspace(-3,1,30),0.683);
h1=loglog(xm,ym,'r','linewidth',3);hold on;
plot(xm,yl,'r--','linewidth',3);
l=legend([h1,h2],'Data','Mock');set(l,'location','southeast','interpreter','latex');
xlabel('$R_{50}[\rm{Mpc}/h]$');
ylabel('$L_{grp}$[L$_\odot/h^2$]');
plot([5e-3,5e-3],[1e9,1e13],'k--');
xlim([1e-3,10]);ylim([1e9,1e13]);
% plot(grpmock.Rad50(grpmock.Mult>2),grpmock.Luminosity(grpmock.Mult>2),'g.');
% plot(grp.Rad50(grp.Mult>2),grp.Luminosity(grp.Mult>2),'r.');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/L-R.eps');
%%
myfigure;
for vol=1:9
    load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
f=grpmock.Mult>2;
[xm,ym,yl]=skeleton(grpmock.VelDisp(f),grpmock.Luminosity(f),logspace(1.5,3.5,25),0.683);
h2=loglog(xm,ym,'g','linewidth',1);hold on;
end
plot(xm,yl,'g--','linewidth',1);
f=grp.Mult>2;
[xm,ym,yl]=skeleton(grp.VelDisp(f),grp.Luminosity(f),logspace(1.5,3.5,25),0.683);
h1=loglog(xm,ym,'r','linewidth',3);hold on;
plot(xm,yl,'r--','linewidth',3);
l=legend([h1,h2],'Data','Mock');set(l,'location','southeast','interpreter','latex');
xlabel('$\sigma_v$[km/s]');
ylabel('$L_{grp}$[L$_\odot/h^2$]');
ylim([1e9,1e13]);xlim([10,1e4]);
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/L-Sigma.eps');
%%
f=grp.Mult>2;
[xm,ym,yl]=skeleton(grp.Rad50(f),grp.VolMult(f),logspace(-4,1,30),0.683);
figure;
loglog(xm,ym,'r');hold on;
plot(xm,yl,'r--');
f=grpmock.Mult>2;
[xm,ym,yl]=skeleton(grpmock.Rad50(f),grpmock.VolMult(f),logspace(-4,1,30),0.683);
loglog(xm,ym,'g');
plot(xm,yl,'g--');
%%
f=grp.Mult>2;
[xm,ym,yl]=skeleton(grp.VelDisp(f),grp.VolMult(f),logspace(1.5,3.5,30),0.683);
figure;
loglog(xm,ym,'r');hold on;
plot(xm,yl,'r--');
f=grpmock.Mult>2;
[xm,ym,yl]=skeleton(grpmock.VelDisp(f),grpmock.VolMult(f),logspace(1.5,3.5,30),0.683);
loglog(xm,ym,'g');
plot(xm,yl,'g--');
%%
f=grp.Mult>2;
[xm,ym,yl]=skeleton(grp.Luminosity(f),grp.VolMult(f),logspace(8,12,30),0.683);
figure;
loglog(xm,ym,'r');hold on;
plot(xm,yl,'r--');
f=grpmock.Mult>2;
[xm,ym,yl]=skeleton(grpmock.Luminosity(f),grpmock.VolMult(f),logspace(8,12,30),0.683);
loglog(xm,ym,'g');
plot(xm,yl,'g--');
%%
f=grp.Mult>2;
[xm,ym,yl]=skeleton(grp.VolMult(f),grp.Luminosity(f),logspace(1,5,30),0.683);
figure;
loglog(xm,ym,'r');hold on;
plot(xm,yl,'r--');
f=grpmock.Mult>2;
[xm,ym,yl]=skeleton(grpmock.VolMult(f),grpmock.Luminosity(f),logspace(1,5,30),0.683);
loglog(xm,ym,'g');
plot(xm,yl,'g--');
%%
f=grp.Mult>2;
[xm,ym,yl]=skeleton(grp.Rad50(f),grp.VelDisp(f),logspace(-4,1,30),0.683);
figure;
loglog(xm,ym,'r');hold on;
plot(xm,yl(:,1),'r--');
plot(xm,yl(:,2),'r--');
f=grpmock.Mult>2;
[xm,ym,yl]=skeleton(grpmock.Rad50(f),grpmock.VelDisp(f),logspace(-4,1,30),0.683);
loglog(xm,ym,'g');
plot(xm,yl(:,1),'g--');
plot(xm,yl(:,2),'g--');
%%
[xm,ym,yl]=skeleton(grp.Rad50,grp.IterCenSM,logspace(-4,1,30),0.683);
figure;
loglog(xm,ym,'r');hold on;
plot(xm,yl(:,1),'r--');
plot(xm,yl(:,2),'r--');
[xm,ym,yl]=skeleton(grpmock.Rad50,grpmock.MstarIter,logspace(-4,1,30),0.683);
loglog(xm,ym,'g');
plot(xm,yl(:,1),'g--');
plot(xm,yl(:,2),'g--');
%%
[xm,ym,yl]=skeleton(grp.IterCenSM,grp.Rad50,logspace(8,12,30),0.683);
figure;
loglog(xm,ym,'r');hold on;
plot(xm,yl(:,1),'r--');
plot(xm,yl(:,2),'r--');
[xm,ym,yl]=skeleton(grpmock.MstarIter,grpmock.Rad50,logspace(8,12,30),0.683);
loglog(xm,ym,'g');
plot(xm,yl(:,1),'g--');
plot(xm,yl(:,2),'g--');
%%
[xm,ym,yl]=skeleton(grp.Luminosity,grp.Rad50,logspace(8,12,30),0.683);
figure;
loglog(xm,ym,'r');hold on;
plot(xm,yl(:,1),'r--');
plot(xm,yl(:,2),'r--');
[xm,ym,yl]=skeleton(grpmock.Luminosity,grpmock.Rad50,logspace(8,12,30),0.683);
loglog(xm,ym,'g');
plot(xm,yl(:,1),'g--');
plot(xm,yl(:,2),'g--');
%%
f=grp.Mult>2;
[xm,ym,yl]=skeleton(grp.VelDisp(f),grp.Rad50(f),logspace(1,4,30),0.683);
figure;
loglog(xm,ym,'r');hold on;
plot(xm,yl,'r--');
f=grpmock.Mult>2;
[xm,ym,yl]=skeleton(grpmock.VelDisp(f),grpmock.Rad50(f),logspace(1,4,30),0.683);
loglog(xm,ym,'g');
plot(xm,yl,'g--');
%%
f=grp.d3area>0;
[xm,ym,yl]=skeleton(grp.Luminosity(f),grp.d3area(f),logspace(8,12,30),0.683);
figure;
loglog(xm,ym,'r');hold on;
plot(xm,yl(:,1),'r--');
plot(xm,yl(:,2),'r--');
f=grpmock.d3area>0;
[xm,ym,yl]=skeleton(grpmock.Luminosity(f),grpmock.d3area(f),logspace(8,12,30),0.683);
loglog(xm,ym,'g');
%%
f=grp.d3area>0;
[xm,ym,yl]=skeleton(grp.d3area(f),grp.Luminosity(f),logspace(-3,2,30),0.683);
figure;
loglog(xm,ym,'r');hold on;
plot(xm,yl(:,1),'r--');
plot(xm,yl(:,2),'r--');
f=grpmock.d3area>0;
[xm,ym,yl]=skeleton(grpmock.d3area(f),grpmock.Luminosity(f),logspace(-3,2,30),0.683);
loglog(xm,ym,'g');
plot(xm,yl(:,1),'g--');
plot(xm,yl(:,2),'g--');
%%
f=grp.d3vol>0;
[xm,ym,yl]=skeleton(grp.Luminosity(f),grp.d3vol(f),logspace(8,12,30),0.683);
figure;
loglog(xm,ym,'r');hold on;
plot(xm,yl(:,1),'r--');
plot(xm,yl(:,2),'r--');
f=grpmock.d3vol>0;
[xm,ym,yl]=skeleton(grpmock.Luminosity(f),grpmock.d3vol(f),logspace(8,12,30),0.683);
loglog(xm,ym,'g');
plot(xm,yl(:,1),'g--');
plot(xm,yl(:,2),'g--');
%%
f=grp.d3vol>0;
[xm,ym,yl]=skeleton(grp.d3vol(f),grp.Luminosity(f),logspace(-3,2,30),0.683);
figure;
loglog(xm,ym,'r');hold on;
plot(xm,yl(:,1),'r--');
plot(xm,yl(:,2),'r--');
f=grpmock.d3vol>0;
[xm,ym,yl]=skeleton(grpmock.d3vol(f),grpmock.Luminosity(f),logspace(-3,2,30),0.683);
loglog(xm,ym,'g');
plot(xm,yl(:,1),'g--');
plot(xm,yl(:,2),'g--');
%%
f=grp.d2radec>0;
[xm,ym,yl]=skeleton(grp.d2radec(f),grp.Luminosity(f),logspace(-3,2,30),0.683);
figure;
loglog(xm,ym,'r');hold on;
plot(xm,yl(:,1),'r--');
plot(xm,yl(:,2),'r--');
f=grpmock.d2radec>0;
[xm,ym,yl]=skeleton(grpmock.d2radec(f),grpmock.Luminosity(f),logspace(-3,2,30),0.683);
loglog(xm,ym,'g');
plot(xm,yl(:,1),'g--');
plot(xm,yl(:,2),'g--');
%%
f=grp.Mult>2;
[xm,ym,yl]=skeleton(grp.VolMult(f),grp.VelDisp(f),logspace(1.3,4.3,30),0.683);
figure;
loglog(xm,ym,'r');hold on;
plot(xm,yl(:,1),'r--');
plot(xm,yl(:,2),'r--');
f=grpmock.Mult>2;
[xm,ym,yl]=skeleton(grpmock.VolMult(f),grpmock.VelDisp(f),logspace(1.3,4.3,30),0.683);
loglog(xm,ym,'g');
plot(xm,yl(:,1),'g--');
plot(xm,yl(:,2),'g--');

%%
f=grp.Mult>2;
[xm,ym,yl]=skeleton(grp.Mult(f),grp.VelDisp(f),logspace(0,3,30),0.683);
figure;
loglog(xm,ym,'r');hold on;
plot(xm,yl(:,1),'r--');
plot(xm,yl(:,2),'r--');
f=grpmock.Mult>2;
[xm,ym,yl]=skeleton(grpmock.Mult(f),grpmock.VelDisp(f),logspace(0,3,30),0.683);
loglog(xm,ym,'g');
plot(xm,yl(:,1),'g--');
plot(xm,yl(:,2),'g--');
%%
f=grp.Mult>2;
[xm,ym,yl]=skeleton(grp.Luminosity(f),grp.Mult(f),logspace(8,12,30),0.683);
figure;
loglog(xm,ym,'r');hold on;
plot(xm,yl(:,1),'r--');
plot(xm,yl(:,2),'r--');
f=grpmock.Mult>2;
[xm,ym,yl]=skeleton(grpmock.Luminosity(f),grpmock.Mult(f),logspace(8,12,30),0.683);
loglog(xm,ym,'g');
plot(xm,yl(:,1),'g--');
plot(xm,yl(:,2),'g--');
%%
nbin=30;
f=grp.Mult>2;
figure;
[xx,yy,n,s]=densitygrid(log10(grp.Rad50(f)),log10(grp.VolMult(f)),[nbin,nbin],[-3,2],[1.3,4.3]);
[~,h2]=contour(xx,yy,log10((n+1)),1:0.4:3,'r');hold on;
f=grpmock.Mult>2;
[xx,yy,n,s]=densitygrid(log10(grpmock.Rad50(f)),log10(grpmock.VolMult(f)),[nbin,nbin],[-3,2],[1.3,4.3]);
[~,h2]=contour(xx,yy,log10((n+1)),1:0.4:3,'b');hold on;
%%
nbin=30;
f=grp.Mult>2;
figure;
[xx,yy,n,s]=densitygrid(log10(grp.Luminosity(f)),log10(grp.VolMult(f)),[nbin,nbin],[8,12],[1.3,4.3]);
[~,h2]=contour(xx,yy,log10((n+1)),1:0.4:3,'r');hold on;
f=grpmock.Mult>2;
[xx,yy,n,s]=densitygrid(log10(grpmock.Luminosity(f)),log10(grpmock.VolMult(f)),[nbin,nbin],[8,12],[1.3,4.3]);
[~,h2]=contour(xx,yy,log10((n+1)),1:0.4:3,'b');hold on;
%%
nbin=50;
f=grp.Mult>2;
figure;
[xx,yy,n,s]=densitygrid(log10(grp.Luminosity(f)),log10(grp.Rad50(f)),[nbin,nbin],[8,13],[-3,2]);
[~,h2]=contour(yy,xx,log10((n+1)),6,'r');hold on;
f=grpmock.Mult>2;
[xx,yy,n,s]=densitygrid(log10(grpmock.Luminosity(f)),log10(grpmock.Rad50(f)),[nbin,nbin],[8,13],[-3,2]);
[~,h2]=contour(yy,xx,log10((n+1)),6,'b');hold on;