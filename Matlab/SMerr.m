nbin=50;
smth=@(n) conv2(n,[0.05,0.1,0.05;0.1,0.4,0.1;0.05,0.1,0.05;],'same');
% smth=@(n) conv2(n,[0.08,0.15,0.08;0.15,0.23,0.15;0.08,0.15,0.08;],'same');
percents=[0.3,0.6,0.9];
myfigure;
f=gal.CentralSampleIter>0;
x=log10(gal.SMsps(f));y=(gal.SMerr(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[9,12],[0.08,0.2]);
levels=percentile_to_density(n,percents);
[~,h0]=contour(xx,yy,smth(n),levels,'r');hold on;%colormap(cmap);
[xm,ym,yl]=skeleton(x,y,9:0.2:11.8);
plot(xm,ym,'k-');plot(xm,yl,'k--');
xlabel('$M_\star$[M$_\odot/h^2$]');
ylabel('$\sigma_{\log M_\star}$');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/extra/SMerr.eps');
%%
f=gal.CentralSampleIter>0;
x=log10(gal.SMsps(f));y=(gal.SMerr(f));
[xm,ym,yl]=skeleton(x,y,9:0.3:12);
plot(xm,ym,'b-');
hold on;
plot(xm,yl,'b--');
%%

nbin=50;
smth=@(n) conv2(n,[0.05,0.1,0.05;0.1,0.4,0.1;0.05,0.1,0.05;],'same');
% smth=@(n) conv2(n,[0.08,0.15,0.08;0.15,0.23,0.15;0.08,0.15,0.08;],'same');
percents=[0.3,0.6,0.9];
% myfigure;
f=gal.CentralSampleIter>0&gal.SSFR>10^-1.5&gal.FlagSFR>0;
x=log10(gal.SMsps(f));y=(gal.SMerr(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[9,12],[0.08,0.2]);
levels=percentile_to_density(n,percents);
[~,h0]=contour(xx,yy,smth(n),levels,'g');hold on;%colormap(cmap);

[xm,ym,yl]=skeleton(x,y,9:0.2:11.8);
plot(xm,ym,'g--');
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/SMerr.eps');