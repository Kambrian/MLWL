figure;
hold on;
[h,r1,s1,es1]=plotgama_prof('H1','r.');
hold on;
[h,r2,s2,es2]=plotgama_prof('H2','bo');
[h,r3,s3,es3]=plotgama_prof('H3','ms');
%%
figure;[h,r1c,s1c,es1c]=plotgama_prof('H1.CC2GAMA','r.');
hold on;
[h,r2c,s2c,es2c]=plotgama_prof('H2.CC2GAMA','go');
[h,r3c,s3c,es3c]=plotgama_prof('H3.CC2GAMA','bs');
%%
myfigure;
% [h,r,s,es]=plotgama_prof('H1.LONG','r.');
% hold on;
% [h,r,s,es]=plotgama_prof('H2.LONG','go');
% [h,r,s,es]=plotgama_prof('H3.LONG','r.');


name='H3.LONG';
file=[datadir,'randprof_',name,'.txt'];
rand=importdata(file,'\t',1);
file=[datadir,'srotateprof_',name,'.txt'];
srotate=importdata(file,'\t',1);
file=[datadir,'shearprof_',name,'.txt'];
shear=importdata(file,'\t',1);
file=[datadir,'shearcov_',name,'.txt'];
cov=load(file);
file=[datadir,'srotatecov_',name,'.txt'];
srotatecov=load(file);
p=shear.data(:,1);
s=shear.data(:,2)./p;
es=sqrt(diag(cov))./p;
errorbar(shear.data(:,5),signlog(s),signlog(s)-signlog(s-es),signlog(s+es)-signlog(s),'r.')
set(gca,'xscale','log');
hold on;
p=srotate.data(:,1);
s=srotate.data(:,2)./p;
es=srotate.data(:,6)./p./sqrt(srotate.data(:,end));
errorbar(shear.data(:,5),signlog(s),signlog(s)-signlog(s-es),signlog(s+es)-signlog(s),'gs')
p=rand.data(:,1);
s=rand.data(:,2)./p;
es=rand.data(:,6)./p./sqrt(rand.data(:,end));
errorbar(shear.data(:,5),signlog(s),signlog(s)-signlog(s-es),signlog(s+es)-signlog(s),'bo')
legend('GAMA','ShearRotate','RandLens');

xlabel('r/(Mpc/h)');
ylabel('$slog(\Delta\Sigma)$');
print('-depsc','wlSig_H3_LONG.eps');
%%
myfigure;
% [h,r,s,es]=plotgama_prof('H1.LONG','r.');
% hold on;
% [h,r,s,es]=plotgama_prof('H2.LONG','go');
[h,r,s,es]=plotgama_prof('H3.LONG','ro');
xlabel('r/(Mpc/h)');
ylabel('$slog(\Delta\Sigma)$');
print('-depsc','wlSig_H3_LONG_subtract.eps');
%%
figure;
[h,r,w,ew]=plotgama_corr('H1','r.');
hold on;
[h,r,w,ew]=plotgama_corr('H2','go');
[h,r,w,ew]=plotgama_corr('H3','bs');
figure;[h,r,w,ew]=plotgama_corr('H1.CC2GAMA','r.');
hold on;
[h,r,w,ew]=plotgama_corr('H2.CC2GAMA','go');
[h,r,w,ew]=plotgama_corr('H3.CC2GAMA','bs');

%%
myfigure;
figure;
[h,r,w,ew]=plotgama_corr('H1','r.');hold on;
[h,r,w1,ew1]=plotgama_corr('H1.CC2GAMA','go');
close;
% w=w-1;w1=w1-1;
% errorbar(r,w1./w,sqrt(ew1.^2./w1.^2+ew.^2./w.^2).*(w1./w),'r.');
errorbar(r,w1-w,sqrt(ew1.^2+ew.^2),'r.');
hold on;
figure;
[h,r,w,ew]=plotgama_corr('H2','r.');
[h,r,w1,ew1]=plotgama_corr('H2.CC2GAMA','go');
close;
% w=w-1;w1=w1-1;
% errorbar(r,w1./w,sqrt(ew1.^2./w1.^2+ew.^2./w.^2).*(w1./w),'go');
errorbar(r,w1-w,sqrt(ew1.^2+ew.^2),'go');
figure;
[h,r,w,ew]=plotgama_corr('H3','r.');
[h,r,w1,ew1]=plotgama_corr('H3.CC2GAMA','go');
close;
% w=w-1;w1=w1-1;
% errorbar(r,w1./w,sqrt(ew1.^2./w1.^2+ew.^2./w.^2).*(w1./w),'bs');
errorbar(r,w1-w,sqrt(ew1.^2+ew.^2),'bs');
plot([1e-2,10],[0,0],'k');
set(gca,'xscale','log');
ylim([-0.2,0.1]);
xlabel('r/(Mpc/h)');
ylabel('$B_{CC2}\frac{\ }{\ }B_{ZEBRA}$');
l=legend('M_{Low}','M_{Mid}','M_{Hig}'); set(l,'location','southeast');
print('-depsc','/home/kam/Projects/Lensing/output/boost_cmp_CC2ZEBRA.eps');

%% shearprof to average photoz prof
myfigure;
plotgama_prof('H3.LONG','ro--',2)
legend('Observation','Photo-z Average')
print('-depsc','photoz_av_H3_LONG.eps');