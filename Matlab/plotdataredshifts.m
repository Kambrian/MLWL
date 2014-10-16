%load_G3Cv4
load G3Cv4up8/G3Cv4up8.mat
zgrp=grp.Zfof;
zgal=gal.Zspec;

f=grp.Mult>2;

load /work/Projects/Lensing/data/shearmap_ZEBRA.mat
zsrc=[e1(:,6);e2(:,6);e3(:,6)];
%%
myfigure;
[xs,ys,~,h]=linhist(zsrc,20,'stairs','k-');
% set(h,'displayname',['source']);
hold on;
[xg2,yg2,~,h]=linhist(zgrp(f),10,'stairs','k--');
% set(h,'displayname',['groups']);
set(gca,'yscale','log');
xmax=max(zgrp(f))+2*eps;
plot([xmax,xmax],[0.1,yg2(end)],'k--');
ylim([1,1e6]);
xlim([0,1.15]);
l=legend('source','groups'); set(l,'interpreter','latex','fontsize',25,'position',[0.74,0.82,0.1,0.1]);
xlabel('Redshift');ylabel('Counts');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/RedshiftFun.eps');
%%
myfigure;
[xs,ys,~,h]=linhist(zsrc,20,'stairs','r-');
set(h,'displayname',['source,',num2str(numel(zsrc))]);
hold on;
[xg,yg,~,h]=linhist(zgrp,10,'stairs','g-');
set(h,'displayname',['groups,',num2str(numel(zgrp))]);
% plot(xs,ys,'r-');
hold on;
% plot(xg,yg,'g-');
[xg2,yg2,~,h]=linhist(zgrp(f),10,'stairs','g--');
set(h,'displayname',['groups(Mult>2),',num2str(numel(zgrp(f)))]);
% plot(xg2,yg2,'g--');
[xgal,ygal,~,h]=linhist(zgal,10,'stairs','k-');
set(h,'displayname',['GAMA gals,',num2str(numel(zgal))]);
set(gca,'yscale','log');
legend('show');
xlabel('z');ylabel('dN');
% print('-depsc','DATAredshifts.eps');