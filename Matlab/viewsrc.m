data=loadShearTBL(1);
er=data(:,28:30);
ei=data(:,34:36);
ea=data(:,40:42);
r=data(:,16);
%% PSF-size map
ncolor=256;
myfigure;
c=colormap(jet(ncolor));
ra=data(:,1);
dec=data(:,2);
y=data(:,33);
[tmp,ic]=histc(y,linspace(min(y)-eps,max(y)+eps,ncolor));
%      [tmp,ic]=histc(y,linspace(-5,5,ncolor));    ic(data>5)=ncolor;ic(data<-5)=1;
for j=1:10:numel(ra)
    plot(ra(j),dec(j),'.','color',c(ic(j),:));
    hold on;
end
axis tight;
print('-depsc','PSF_map_1.eps');
%%
x=10:0.1:21.8;
figure;[xm,ym,dyx,h,xr]=linhist(r,x,'stairsnorm');
y=cumsum(ym);
figure;plot(xr(105:end),log10(y(105:end)),'.')
x=16:0.02:21.8;
[xm,ym,dyx,h,xr]=linhist(r,x);
figure;plot(xm,log10(dyx),'.');
%%
err=ea(:,3)*2; % measurement noise of e, single component? 
signal=ea(:,1);
signal2=ea(:,2);
x=16:0.2:21.8;
[n,bin]=histc(r,x);
x=x(1:end-1);
xm=x;ya=zeros(numel(x),1);yr=ya;
sbrd=zeros(numel(x),2);
smed=ya;svar=ya;emed=ya;
evar=emed;
alpha=(1-0.683)/2;
myfigure;
colors=colormap(jet(numel(x)));
for i=1:numel(x)
    if n(i)
    xm(i)=median(r(bin==i));
    yr(i)=median(er(bin==i,3));
    tmp=sort(signal(bin==i));
    sbrd(i,:)=[tmp(ceil(alpha*n(i))),tmp(ceil((1-alpha)*n(i)))];
    smed(i)=median(signal(bin==i));
    svar(i)=sqrt(std(signal(bin==i)).^2+std(signal2(bin==i)).^2)/sqrt(2);
    emed(i)=median(err(bin==i));
    evar(i)=sqrt(mean(err(bin==i).^2));
    [tmp1,tmp2]=linhist(signal(bin==i),-2:0.2:2);
    plot(tmp1,tmp2/sum(tmp2),'--','color',colors(i,:));hold on;
    end
end
xlabel('ellipticity');ylabel('Probability');colorbar;
title('bin in r model mag');
print('-depsc','Ellip_distr_magbin.eps');
figure;plot(r,signal,'.');
hold on;
plot(xm,smed,'r-');
plot(xm,smed+svar,'g-',xm,smed-svar,'g-');
plot(xm,sbrd(:,1),'k--',xm,sbrd(:,2),'k--');

myfigure;
% plot(r,ea(:,3)*2,'.');hold on;
plot(xm,svar,'r-');hold on;
a=sqrt(svar.^2-evar.^2);
plot(xm,a,'k-');
plot(xm,evar,'g-');
xlabel('r model magnitude');ylabel('single component $\sigma_e$');
eav=sqrt(sum((a.^2.*n(1:end-1)))./sum(n(1:end-1))); %intrinsic noise
plot(xm,repmat(eav,size(xm)),'--');
legend('Total variance','Intrinsic','Measurement','intrinsic average');
% print('-depsc','/home/kam/Projects/Lensing/output/shape_noise.eps');

myfigure;
[xx,yy,n,s]=densitygrid(r(r>16),signal(r>16),[20,20]);
mesh(xx,yy,n)
% contourf(xx,yy,n);
xlabel('$r$');
ylabel('$e$');
% print('-depsc','Ellip_surf_magbin.eps');
%%
err=ea(:,3)*2; % measurement noise of e, single component? 
signal=ea(:,1);
signal2=ea(:,2);
x=0:0.05:0.6;
[n,bin]=histc(err,x);
x=x(1:end-1);
xm=x;ya=zeros(numel(x),1);yr=ya;
sbrd=zeros(numel(x),2);
smed=ya;svar=ya;emed=ya;
evar=emed;
alpha=(1-0.683)/2;
myfigure;
colors=colormap(jet(numel(x)));
for i=1:numel(x)
    if n(i)
    xm(i)=median(err(bin==i));
    yr(i)=median(er(bin==i,3));
    tmp=sort(signal(bin==i));
    sbrd(i,:)=[tmp(ceil(alpha*n(i))),tmp(ceil((1-alpha)*n(i)))];
    smed(i)=median(signal(bin==i));
    svar(i)=sqrt(std(signal(bin==i)).^2+std(signal2(bin==i)).^2)/sqrt(2);
    emed(i)=median(err(bin==i));
    evar(i)=sqrt(mean(err(bin==i).^2));
    [tmp1,tmp2]=linhist(signal(bin==i),-2:0.2:2);
    plot(tmp1,tmp2/sum(tmp2),'--','color',colors(i,:));hold on;
    end
end
xlabel('ellipticity');ylabel('Probability');colorbar;
title('bin in measurement error');
% print('-depsc','Ellip_distr_errbin.eps');
figure;plot(err,signal,'.');
hold on;
plot(xm,smed,'r-');
plot(xm,smed+svar,'g-',xm,smed-svar,'g-');
plot(xm,sbrd(:,1),'k--',xm,sbrd(:,2),'k--');

myfigure;
% plot(err,ea(:,3),'.');hold on;
plot(xm,svar,'r-');hold on;
a=sqrt(svar.^2-evar.^2);
plot(xm,a,'k-');
plot(xm,evar,'g-');
xlabel('measurement error');ylabel('single component $\sigma_e$');
eav=sqrt(sum((a.^2.*n(1:end-1)))./sum(n(1:end-1))); %intrinsic noise
plot(xm,repmat(eav,size(xm)),'--');
legend('Total variance','Intrinsic','Measurement','intrinsic average');
% print('-depsc','/home/kam/Projects/Lensing/output/shape_noise.eps');

myfigure;
[xx,yy,n,s]=densitygrid(err(err<0.6),signal(err<0.6),[20,20]);
mesh(xx,yy,n)
% contourf(xx,yy,n);
set(gca,'xlim',[0,0.6]);
xlabel('$measurement error$');
ylabel('$ellipticity$');
% print('-depsc','Ellip_surf_errbin.eps');
%% total noise versus measurement noise
nbin=50;
err=ea(:,3);
signal=ea(:,1);

emin=min(err);
emax=max(err);
% ebin=logspace(log10(emin),log10(emax),nbin+1);
ebin=linspace((emin),(emax),nbin+1);
[n,bin]=histc(err,ebin);
smed=zeros(nbin,1);svar=smed;emed=smed;
sbrd=zeros(nbin,2);
alpha=(1-0.683)/2;
for i=1:nbin
    if n(i)
    tmp=sort(signal(bin==i));
    sbrd(i,:)=[tmp(ceil(alpha*n(i))),tmp(ceil((1-alpha)*n(i)))];
    smed(i)=median(signal(bin==i));
    svar(i)=std(signal(bin==i));
    emed(i)=median(err(bin==i));
    end
end
figure;plot(err,signal,'.');
hold on;
plot(emed,smed,'r-');
plot(emed,smed+svar,'g-',emed,smed-svar,'g-');
plot(emed,sbrd(:,1),'k--',emed,sbrd(:,2),'k--');
figure(2); 
%hold on;
plot(emed,svar,'ro');
axis equal
%% err vs. z
z=data(:,53);
err=ea(:,3)*2; % measurement noise of e, single component? 
signal=ea(:,1); % average ellipticity of r and i band, first component
signal2=ea(:,2);
x=0:0.1:1;
[n,bin]=histc(z,x);
x=x(1:end-1);
xm=x;ya=zeros(numel(x),1);yr=ya;
sbrd=zeros(numel(x),2);
smed=ya;svar=ya;emed=ya;
evar=emed;
alpha=(1-0.683)/2;
myfigure;
colors=colormap(jet(numel(x)));
for i=1:numel(x)
    if n(i)
    xm(i)=median(z(bin==i));
    yr(i)=median(er(bin==i,3));
    tmp=sort(signal(bin==i));
    sbrd(i,:)=[tmp(ceil(alpha*n(i))),tmp(ceil((1-alpha)*n(i)))];
    smed(i)=median(signal(bin==i));
    svar(i)=sqrt(std(signal(bin==i)).^2+std(signal2(bin==i)).^2)/sqrt(2);
    emed(i)=median(err(bin==i));
    evar(i)=sqrt(mean(err(bin==i).^2));
    [tmp1,tmp2]=linhist(signal(bin==i),-2:0.2:2);
    plot(tmp1,tmp2/sum(tmp2),'--','color',colors(i,:));hold on;
    end
end
xlabel('ellipticity');ylabel('Probability');colorbar;
title('bin in redshift');
% print('-depsc','Ellip_distr_magbin.eps');
figure;plot(z,signal,'.');
hold on;
plot(xm,smed,'r-');
plot(xm,smed+svar,'g-',xm,smed-svar,'g-');
plot(xm,sbrd(:,1),'k--',xm,sbrd(:,2),'k--');

myfigure;
% plot(r,ea(:,3)*2,'.');hold on;
plot(xm,svar,'r-');hold on;
a=sqrt(svar.^2-evar.^2);
plot(xm,a,'k-');
plot(xm,evar,'g-');
xlabel('z');ylabel('single component $\sigma_e$');
eav=sqrt(sum((a.^2.*n(1:end-1)))./sum(n(1:end-1))); %intrinsic noise
plot(xm,repmat(eav,size(xm)),'--');
legend('Total variance','Intrinsic','Measurement','intrinsic average');
% print('-depsc','/home/kam/Projects/Lensing/output/shape_noise.eps');

myfigure;
[xx,yy,n,s]=densitygrid(z(z<1),signal(z<1),[20,20]);
mesh(xx,yy,n)
% contourf(xx,yy,n);
xlabel('$z$');
ylabel('$e$');
% print('-depsc','Ellip_surf_magbin.eps');
%%
z=data(:,53);

[xm,ym,dyx]=linhist(z,30);
figure;plot(xm,dyx/sum(ym),'--.');
hold on;
plot(xm,source_zdistr(xm))


%%
x=-2:0.2:2;
figure;linhist(er(:,1),x,'stairs','r');
hold on;
linhist(ei(:,1),x,'stairs','g');
linhist(ea(:,1),x,'stairs','k');

figure;linhist(er(:,2),x,'stairs','r');
hold on;
linhist(ei(:,2),x,'stairs','g');
linhist(ea(:,2),x,'stairs','k');

x=0:0.1:1;
figure;linhist(er(:,3),x,'stairsnorm','r');
hold on;
linhist(ei(:,3),x,'stairsnorm','g');
linhist(ea(:,3),x,'stairsnorm','k');

figure;plot((ei(:,1)+er(:,1))/2,ea(:,1),'.')
