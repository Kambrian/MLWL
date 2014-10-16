function compare_masscontour(grpmock,ma,mb,c)

if nargin<4
    c='k';
end
nbin=50;
mlim=[11,15];
alphaplot=[5,95];
MultCut=5;
smth=@(n) conv2(n,[0.05,0.1,0.05;0.1,0.4,0.1;0.05,0.1,0.05;],'same');
percents=[0.3,0.6,0.9];
xbin=linspace(11,15,15);
f=grpmock.Mult>2;

x=log10(grpmock.(ma));
y=log10(grpmock.(mb));

[xx,yy,n]=densitygrid(x(f),y(f),[nbin,nbin],mlim,mlim);
h0=contour(xx,yy,smth(n),percentile_to_density(n,percents),c);hold on;
% ff=grpmock.Zfof<0.2&grpmock.Mult<5&grpmock.Mult>2;
% [xx,yy,n]=densitygrid(x(ff),y(ff),[nbin,nbin],mlim,mlim);
% contour(xx,yy,smth(n),percentile_to_density(n,percents),'r');hold on;
% ff=grpmock.Zfof>0.2&grpmock.Mult>5;
% [xx,yy,n]=densitygrid(x(ff),y(ff),[nbin,nbin],mlim,mlim);
% contour(xx,yy,smth(n),percentile_to_density(n,percents),'g');hold on;
ff=grpmock.Mult>2;
[xmed,ymed,yl]=skeleton(x(ff),y(ff),xbin);
h1=plot(xmed,ymed,'r-');
plot(xmed,yl,'r-');
% 
% ff=grpmock.Zfof<0.2&grpmock.Mult>2;
% xbin=prctile(x(ff),alphaplot);
% xbin=linspace(xbin(1),xbin(2),10);
% [xmed,ymed]=skeleton(x(ff),y(ff),xbin);
% h2=plot(xmed,ymed,'r--');
% ff=grpmock.Zfof>0.2&grpmock.Mult>2;
% xbin=prctile(x(ff),alphaplot);
% xbin=linspace(xbin(1),xbin(2),10);
% [xmed,ymed]=skeleton(x(ff),y(ff),xbin);
% h3=plot(xmed,ymed,'g--');
% ff=grpmock.Mult>2&grpmock.Mult<MultCut;
% xbin=prctile(x(ff),alphaplot);
% xbin=linspace(xbin(1),xbin(2),10);
% [xmed,ymed]=skeleton(x(ff),y(ff),xbin);
% h4=plot(xmed,ymed,'r-');
% ff=grpmock.Mult>=MultCut;
% xbin=prctile(x(ff),alphaplot);
% xbin=linspace(xbin(1),xbin(2),10);
% [xmed,ymed]=skeleton(x(ff),y(ff),xbin);
% h5=plot(xmed,ymed,'g-');
ylim(mlim);xlim(mlim);plot(mlim,mlim,'k:');
xlabel(ma);
ylabel(mb);
% text(11.2,14.5,ma,'fontsize',15);