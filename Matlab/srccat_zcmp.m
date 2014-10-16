macro=default_WLparam();
macro.ZGAMA=-2;
GAMA_WL_init;
e=[e{1};e{2};e{3}];
zcc2=e(:,6);
err_cc2=sqrt(((e(:,6)-e(:,7)).^2+(e(:,8)-e(:,6)).^2)/2);
macro.ZGAMA=-1;
GAMA_WL_init;
e=[e{1};e{2};e{3}];
zzeb=e(:,6);
err_zeb=sqrt(((e(:,6)-e(:,7)).^2+(e(:,8)-e(:,6)).^2)/2);
macro=default_WLparam();
macro.ZGAMA=2;
GAMA_WL_init;
e=[e{1};e{2};e{3}];
zcc2_gama=e(:,6);
macro.ZGAMA=1;
GAMA_WL_init;
e=[e{1};e{2};e{3}];
zzeb_gama=e(:,6);
macro.ZGAMA=-3;
GAMA_WL_init;
e=[e{1};e{2};e{3}];
zd1=e(:,6);
err_d1=sqrt(((e(:,6)-e(:,7)).^2+(e(:,8)-e(:,6)).^2)/2);
macro.ZGAMA=-4;
GAMA_WL_init;
e=[e{1};e{2};e{3}];
ztmp=e(:,6);
err_tmp=sqrt(((e(:,6)-e(:,7)).^2+(e(:,8)-e(:,6)).^2)/2);
%% matched fraction
cd /home/kam/Projects/Lensing/data
    load src_sdss_match.mat nmatch idmatch zmatch zerrmatch
    zmatch=[zmatch{1};zmatch{2};zmatch{3}];
    nmatch=[nmatch{1};nmatch{2};nmatch{3}];
    matchfraction=[sum(nmatch>0)/numel(nmatch),sum(zmatch>-1)/numel(nmatch)]   %>91%
%%  dN/dz distribution
cd /home/kam/Projects/Lensing/output

myfigure;
[x,y,dyx]=linhist(zcc2,30,'','r--');
hold on;
[x,y2,dyx2]=linhist(zzeb,30,'','g--');
[x,y3,dyx3]=linhist(zcc2_gama,30,'','r-');
[x,y4,dyx4]=linhist(zzeb_gama,30,'','g-');
[x,y5]=linhist(zd1,30,'','b--');
[x,y6]=linhist(ztmp,30,'','c--');
legend('cc2','zebra','cc2+gama','zebra+gama','d1','tmpl');
xlabel('z');ylabel('dN');
print('-depsc','zcmp_srccat.eps');
myfigure;
plot(x,(y-y2),'k.-')
hold on;
plot(x,(y4-y2),'g.-');
plot(x,(y3-y),'r.-');
plot(x,cumsum(y-y2),'-')
l=legend('cc2-zebra','gama-zebra','gama-cc2','cc2-zeb,<z');
set(l,'location','southeast');
xlabel('z');ylabel('$\delta dN$');
print('-depsc','zcmp_srccat_diff.eps');
myfigure;
plot(x,(y-y2)./y3,'k.-')
hold on;
plot(x,(y4-y2)./y3,'g.-');
plot(x,(y3-y)./y3,'r.-');
legend('cc2-zebra','gama-zebra','gama-cc2');
xlabel('z');ylabel('$\delta dN/dN_{gamacc2}$');
% print('-depsc','zcmp_srccat_diff_norm.eps');
%% error dist comparison
f=zcc2>0.4&zcc2<0.5;
g=zzeb>0.4&zzeb<0.5;
figure;
[x,y]=linhist(err_cc2(f)./(1+zcc2(f)),0:0.01:0.5,'stairs','r');  % cc2 peaks at sigma(z)=0.07*(1+z)
% title(num2str(x(y==max(y(x>0.02)))));
hold on;
[x,y]=linhist(err_zeb(g)./(1+zzeb(g)),0:0.01:0.5,'stairs','g');  % cc2 peaks at sigma(z)=0.07*(1+z)
% title(num2str(x(y==max(y(x>0.02)))));
legend('cc2','zebra');
figure;
[x,y]=linhist(err_cc2(f),0:0.01:0.5,'stairs','r');  % cc2 peaks at sigma(z)=0.07*(1+z)
% title(num2str(x(y==max(y(x>0.02)))));
hold on;
[x,y]=linhist(err_zeb(g),0:0.01:0.5,'stairs','g');  % cc2 peaks at sigma(z)=0.07*(1+z)
% title(num2str(x(y==max(y(x>0.02)))));
legend('cc2','zebra');
%% zerr-z distribution
zerr=err_cc2;
z=zcc2;
[xx,yy,n]=densitygrid(z,zerr,[30,30]);
[zmed,zerrmed,~,~,~,~,count]=skeleton(z,zerr,30,0.683);
figure;
contourf(xx,yy,log(n+1));
hold on;
x=0:0.1:1;
plot(x,0.1*x+0.06,'k','linewidth',3);
plot(zmed,zerrmed,'y');
xlabel('z');ylabel('\sigma_z');
title('CC2');
legend('ln(n+1)','\sigma_z=0.1*(z+0.6)','median');
% print('-depsc','zerr_z_CC2.eps');

zerr=err_zeb;
z=zzeb;
[xx,yy,n]=densitygrid(z,zerr,[30,30]);
[zmed,zerrmed,~,~,~,~,count]=skeleton(z,zerr,30,0.683);
figure;
contourf(xx,yy,log(n+1));
hold on;
x=0:0.1:1;
plot(x,0.1*x+0.06,'k','linewidth',3);
plot(zmed,zerrmed,'y');
xlabel('z');ylabel('\sigma_z');
title('ZEBRA');
legend('ln(n+1)','\sigma_z=0.1*(z+0.6)','median');
% print('-depsc','zerr_z_ZEBRA.eps');
