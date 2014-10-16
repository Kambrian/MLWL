clear
clc
load G3Cv4up8/mockcat_1.mat
load G3Cv4up8/G3Cv4up8.mat
%%
ngp=importdata('2dF/redcat.cal0b.genuine.2df_ngp',' ',2);
sgp=importdata('2dF/redcat.cal0b.genuine.2df_sgp',' ',2);
twodf=[ngp.data;sgp.data];
n=twodf(:,8)./twodf(:,5);
z=twodf(:,3);
[xm,ym,yl,~,~,yerr]=skeleton(z,n,0:0.01:0.3,1);
n2df=normrnd(pchip(xm,ym,grpmock.Zfof),pchip(xm,yerr,grpmock.Zfof));
n2df=max(n2df,pchip(xm,yl(:,1),grpmock.Zfof));
n2df=min(n2df,pchip(xm,yl(:,2),grpmock.Zfof));
grpmock.Mult2dF=grpmock.VolMult.*n2df;
grpmock.Mult2dF(grpmock.Zfof>0.3)=0;
grpmock.DynMass2dF=5.0*grpmock.MassProxyRaw;
%% M=A(L)
% data=[3.07611e+10   1.08824e+10  12.27      0.451                   
% 7.08877e+10   1.41219e+10  13.2       0.1728                  
% 1.48753e+11   4.02562e+10  13.49      0.156                   
% 3.6731e+11   9.20875e+10  14.09      0.1561                  
% 8.96239e+11   3.36118e+11  13.98      0.5326                  
% 3.26962e+12   0  14.75      0.5355  ]; %N>2, Luminosity2 (19.8 new calibration, worse than standard).
data=[3.0851e+10   1.08705e+10  12.22      0.5089                  
7.16576e+10   1.39753e+10  13.01      0.2178                  
1.57631e+11   4.07415e+10  13.45      0.1417                  
3.79365e+11   9.89501e+10  13.71      0.2087                  
9.70493e+11   3.28707e+11  14.21      0.2076                  
2.90161e+12   6.07367e+11  14.79      0.3664];
xbrd=[1e10,5e10; 5e10,1e11; 1e11,2.5e11 ;2.5e11,6e11 ;6e11,2e12 ;2e12,5e12];
myfigure;
for vol=1:9
    load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
f=grpmock.Mult>2;
x=logspace(9,13,12);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.Luminosity(f),log10(grpmock.MIter(f)),x,0.68);
h3=plot(xmed,ymed,'g');hold on;
end
plot(xmed,yl(:,1),'g--');
plot(xmed,yl(:,2),'g--');
x=data(:,1);
h1=ploterr(x,data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro','logx');hold on;
x=logspace(10,13);
h2=plot(x,14.23+1.09*log10(x/1e12),'k');
h22=plot(x,14.23-0.07+1.1*log10(x/1e12),'k--');
h222=plot(x,14.23-0.086+(1.08-0.0045)*log10(x/1e12),'b-'); %with 2halo term
% h22=plot(x,14.12+0.97*log10(x/1e12),'k--'); %R<1Rv.
% h22=plot(x,14.20+1.06*log10(x/1e12),'r--'); %excluding z>1 sources.
% h3=plot(x,14.18+1.051*log10(x/1e12),'k--'); %without Xoffset cut
xlabel('$L_{grp}$[L$_\odot/h^2$]');
ylabel('$\log(M_h$[M$_\odot/h$])');
xlim([8e9,6e12]);
set(gca,'xscale','log');
l=legend([h1(1),h2,h22,h222,h3],'Binned','Global','Systematic','2Halo','Mock');set(l,'location','southeast','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/extra/M-L-2Halo.eps');
%% comparison with 2dF
data=[3.0851e+10   1.08705e+10  12.22      0.5089                  
7.16576e+10   1.39753e+10  13.01      0.2178                  
1.57631e+11   4.07415e+10  13.45      0.1417                  
3.79365e+11   9.89501e+10  13.71      0.2087                  
9.70493e+11   3.28707e+11  14.21      0.2076                  
2.90161e+12   6.07367e+11  14.79      0.3664]; 
xbrd=[1e10,5e10; 5e10,1e11; 1e11,2.5e11 ;2.5e11,6e11 ;6e11,2e12 ;2e12,5e12];
myfigure;
% for vol=1
%     load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
f=grpmock.Mult>2;
x=logspace(10,13,12);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.Luminosity(f),log10(grpmock.MIter(f)),x,0.68);
h3=plot(xmed,ymed-log10(xmed),'g');hold on;
% end
plot(xmed,yl-log10([xmed,xmed]),'g--');
ml=[9.2,2.265;9.5,2.2;9.8,2.08;10.1,2.035;10.4,1.95;10.7,2.07;11,2.22;11.3,2.4;11.6,2.48;11.88,2.55];
ml(:,1)=10.^ml(:,1);
plot(ml(:,1),ml(:,2),'-o');hold on;
x=data(:,1);
h1=ploterr(x,data(:,3)-log10(x),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro','logx');hold on;
x=logspace(10,13);
h2=plot(x,14.24+1.09*log10(x/1e12)-log10(x),'k');
xlabel('$L_{grp}$[L$_\odot/h^2$]');
ylabel('$\log(M_h$[M$_\odot/h$])');
xlim([8e9,6e12]);
set(gca,'xscale','log');
l=legend([h1(1),h2,h3],'Binned','Global','Mock');set(l,'location','southeast','interpreter','latex');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/M-L-rat.eps');
%% M/L=A(L)
data=[3.0851e+10   1.08705e+10  13.94      0.3912                  
7.16576e+10   1.39753e+10  14.16      0.2199                  
1.57631e+11   4.07415e+10  14.25      0.1452                  
3.79365e+11   9.89501e+10  14.13      0.2127                  
9.70493e+11   3.28707e+11  14.25      0.2143                  
2.90161e+12   6.07367e+11  14.21      0.3808];
data(:,3)=data(:,3)-12;
xbrd=[1e10,5e10; 5e10,1e11; 1e11,2.5e11 ;2.5e11,6e11 ;6e11,2e12 ;2e12,5e12];
myfigure;
% for vol=1
%     load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
f=grpmock.Mult>2;
x=logspace(9,13,12);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.Luminosity(f),log10(grpmock.MIter(f)./grpmock.Luminosity(f)),x,0.68);
h3=plot(xmed,ymed,'g');hold on;
% end
plot(xmed,yl,'g--');
ml=[9.2,2.265;9.5,2.2;9.8,2.08;10.1,2.035;10.4,1.95;10.7,2.07;11,2.22;11.3,2.4;11.6,2.48;11.88,2.55];
ml(:,1)=10.^ml(:,1);
h4=plot(ml(:,1),ml(:,2),'b-','linewidth',4);hold on;
h1=ploterr(data(:,1),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro','logx');hold on;
xlabel('$L_{grp}$[L$_\odot/h^2$]');
ylabel('$\log(M_h/L_{grp}[h\rm{M}_\odot/L_\odot$])');
set(gca,'xscale','log');
xlim([5e9,6e12]);
ylim([1.2,2.8]);
set(gca,'xtick',[1e10,1e11,1e12]);
l=legend([h1(1),h3,h4],'Binned','Mock','2PIGG');set(l,'location','southeast','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/M-L-rat2.eps');
%% M/L=A(L)
data2dF=[2.65088e+10   1.13231e+10  13.75      0.7392                  
7.11909e+10   1.43482e+10  13.75      0.645                   
1.54719e+11   3.96499e+10  14.49      0.1645                  
3.8212e+11   1.03496e+11  14.45      0.208                   
8.17228e+11   1.99741e+11  14.33      0.341                   
4.33699e+12   0  14.11      0.536 ]; %2dF mult cut
data3=[3.98969e+10   2.38841e+10  13.74      0.489                   
1.54719e+11   3.96499e+10  14.49      0.1645                  
3.8212e+11   1.03496e+11  14.45      0.208                   
8.17228e+11   1.99741e+11  14.33      0.341 ];data3(:,3)=data3(:,3)-12;
data2dF(:,3)=data2dF(:,3)-12;
xbrd=[1e10,5e10; 5e10,1e11; 1e11,2.5e11 ;2.5e11,6e11 ;6e11,2e12 ;2e12,5e12];
myfigure;
f=grpmock.Mult>2;
x=logspace(9,13,12);
% [xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.Luminosity(f),log10(grpmock.MIter(f)./grpmock.Luminosity(f)),x,0.68);
% h3=plot(xmed,ymed,'r');hold on;
f=grpmock.Mult2dF>1;
x=logspace(9,13,20);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.Luminosity(f),log10(grpmock.MIter(f)./grpmock.Luminosity(f)),x,0.68);
h3=plot(xmed,ymed,'g');hold on;
plot(xmed,yl,'g--');
% [xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.Luminosity(f),log10(grpmock.DynMass(f)./grpmock.Luminosity(f)),x,0.68);
% h3=plot(xmed,ymed,'g:');hold on;
ml=[9.2,2.265;9.5,2.2;9.8,2.08;10.1,2.035;10.4,1.95;10.7,2.07;11,2.22;11.3,2.4;11.6,2.48;11.88,2.55];
ml(:,1)=10.^ml(:,1);
h4=plot(ml(:,1),ml(:,2),'b-','linewidth',4);hold on;
x=logspace(10.4,12.3,50);
y=log10(10.^2.4*min(x/10^11.17,1).^1.24);
plot(x,y,'r--');
% h11=ploterr(data2dF(:,1),data2dF(:,3),{xbrd(:,1),xbrd(:,2)},data2dF(:,4),'bd','logx');hold on;
xbrd=[1e10,1e11; 1e11,2.5e11; 2.5e11,6e11; 6e11,2e12];
h1=ploterr(data3(:,1),data3(:,3),{xbrd(:,1),xbrd(:,2)},data3(:,4),'rd','logx');hold on;
x=logspace(10,13);
% h2=plot(x,14.24+1.09*log10(x/1e12)-log10(x),'k');
xlabel('$L_{grp}$[L$_\odot/h^2$]');
ylabel('$\log(M_h/L_{grp}[h\rm{M}_\odot/L_\odot$])');
set(gca,'xscale','log');
xlim([5e9,6e12]);
ylim([1.2,2.8]);
set(gca,'xtick',[1e10,1e11,1e12]);
l=legend([h1(1),h3,h4],'Binned','Mock','2PIGG');set(l,'location','southeast','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/M-L-rat-2PIGGcut.eps');
%%
f=grpmock.Mult2dF>1;
x=logspace(9,13,20);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.Luminosity(f),log10(grpmock.VolMult(f)),x,0.68);
h3=plot(xmed,ymed,'g');hold on;
f=grpmock.Mult>2;
x=logspace(9,13,20);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.Luminosity(f),log10(grpmock.VolMult(f)),x,0.68);
h3=plot(xmed,ymed,'r');hold on;
xscale('log');
yscale('linear');
%%
figure;
f=grpmock.Mult2dF>1&grpmock.Mult>2;
loglog(grpmock.Luminosity(f),grpmock.VolMult(f),'r.');
hold on;
f=grpmock.Mult2dF<2&grpmock.Mult>2;
loglog(grpmock.Luminosity(f),grpmock.VolMult(f),'go');
f=grpmock.Mult2dF>1&grpmock.Mult<=2;
loglog(grpmock.Luminosity(f),grpmock.VolMult(f),'bs');
%%
figure;
f=grpmock.Mult2dF>1&grpmock.Mult>2&grpmock.Luminosity<1e11;
loglog(grpmock.Zfof(f),grpmock.VolMult(f),'r.');
hold on;
f=grpmock.Mult2dF<2&grpmock.Mult>2&grpmock.Luminosity<1e11;
loglog(grpmock.Zfof(f),grpmock.VolMult(f),'g.');
f=grpmock.Mult2dF>1&grpmock.Mult<=2&grpmock.Luminosity<1e11;
loglog(grpmock.Zfof(f),grpmock.VolMult(f),'b.');

f=grpmock.Mult2dF>1&grpmock.Mult>2&grpmock.Luminosity>1e11&grpmock.Luminosity<5e11;
loglog(grpmock.Zfof(f),grpmock.VolMult(f),'rs','markersize',8);
hold on;
f=grpmock.Mult2dF<2&grpmock.Mult>2&grpmock.Luminosity>1e11&grpmock.Luminosity<5e11;
loglog(grpmock.Zfof(f),grpmock.VolMult(f),'gs','markersize',8);
f=grpmock.Mult2dF>1&grpmock.Mult<=2&grpmock.Luminosity>1e11&grpmock.Luminosity<5e11;
loglog(grpmock.Zfof(f),grpmock.VolMult(f),'bs','markersize',8);

f=grpmock.Mult2dF>1&grpmock.Mult>2&grpmock.Luminosity>5e11;
loglog(grpmock.Zfof(f),grpmock.VolMult(f),'rx','markersize',15);
hold on;
f=grpmock.Mult2dF<2&grpmock.Mult>2&grpmock.Luminosity>5e11;
loglog(grpmock.Zfof(f),grpmock.VolMult(f),'gx','markersize',15);
f=grpmock.Mult2dF>1&grpmock.Mult<=2&grpmock.Luminosity>5e11;
loglog(grpmock.Zfof(f),grpmock.VolMult(f),'bx','markersize',15);
%%
figure;
f=grp.Mult2dF>1&grp.Mult>2;
loglog(grp.Luminosity(f),grp.Mult(f),'r.');
hold on;
f=grp.Mult2dF<2&grp.Mult>2;
loglog(grp.Luminosity(f),grp.Mult(f),'go');
f=grp.Mult2dF>1&grp.Mult<=2;
loglog(grp.Luminosity(f),grp.Mult(f),'bs');
%% Understand sample selection
figure;
f=grpmock.Mult>2;
multbin=[2,4,6,10,20,50,100];
massbin=11:0.5:15;
colors=colormap(colorcube(numel(massbin)));
x=grpmock.Mult;y=log10(grpmock.Luminosity);
for i=1:numel(massbin)-1
    ff=f&grpmock.MIter>10.^massbin(i)&grpmock.MIter<=10.^massbin(i+1);
    semilogx(x(ff),y(ff),'o','color',colors(i,:),'markersize',(massbin(i)-8)^1.5);hold on;
end
%%
figure;
f=grpmock.Mult>2;
% multbin=[2,4,6,10,20,50,100];
% massbin=11:0.5:15;
massbin=0:0.1:0.5;
colors=colormap(hot(numel(massbin)));
x=grpmock.Mult;y=log10(grpmock.MIter./grpmock.Luminosity.^2);
for i=1:numel(massbin)-1
%     ff=f&grpmock.MIter>10.^massbin(i)&grpmock.MIter<=10.^massbin(i+1);
ff=f&grpmock.Zfof>massbin(i)&grpmock.Zfof<=massbin(i+1);
    semilogx(x(ff),y(ff),'o','color',colors(i,:),'markersize',massbin(i+1)*10+4);hold on;
end
%%
figure;
f=grpmock.Mult>2;
% multbin=[2,4,6,10,20,50,100];
% massbin=11:0.5:15;
massbin=0:0.1:0.5;
colors=colormap(hot(numel(massbin)));
x=log10(grpmock.Luminosity);y=log10(grpmock.MIter./ffz(grpmock.Zfof));
for i=1:numel(massbin)-1
%     ff=f&grpmock.MIter>10.^massbin(i)&grpmock.MIter<=10.^massbin(i+1);
ff=f&grpmock.Zfof>massbin(i)&grpmock.Zfof<=massbin(i+1);
    plot(x(ff),y(ff),'o','color',colors(i,:),'markersize',massbin(i+1)*10+4);hold on;
end
%%
figure;
f=grpmock.Mult>2;
massbin=[2,4,6,10,20,50,100];
% massbin=11:0.5:15;
% massbin=0:0.1:0.5;
colors=colormap(hot(numel(massbin)));
x=exp(-grpmock.Zfof);
% y=log10(grpmock.MIter./grpmock.Luminosity.^1.8./(grpmock.Zfof).^0);
% y=log10(grpmock.Luminosity);
y=grpmock.Mult./exp(-21.5-14.48*grpmock.Zfof);
for i=1:numel(massbin)-1
%     ff=f&grpmock.MIter>10.^massbin(i)&grpmock.MIter<=10.^massbin(i+1);
ff=f&grpmock.Mult>massbin(i)&grpmock.Mult<=massbin(i+1);
%     plot(log10(ffz(x(ff))),y(ff),'o','color',colors(i,:),'markersize',massbin(i+1)*0.3+4);hold on;
    plot(log10((x(ff))),y(ff),'o','color',colors(i,:),'markersize',massbin(i+1)*0.3+4);hold on;
end
%%
figure;
f=grpmock.Mult>2;
multbin=[2,4,6,10,20,50,100];
massbin=11:0.5:15;
colors=colormap(colorcube(numel(massbin)));
x=grpmock.Mult;y=grpmock.Zfof;
for i=1:numel(massbin)-1
    ff=f&grpmock.MIter>10.^massbin(i)&grpmock.MIter<=10.^massbin(i+1);
    semilogx(x(ff),y(ff),'o','color',colors(i,:),'markersize',(massbin(i)-5)^1.6-15);hold on;
end
%%
figure;
f=grpmock.Mult>2;
multbin=[2,4,6,10,20,50,100];
massbin=9:0.5:13;
colors=colormap(colorcube(numel(massbin)));
x=grpmock.Mult;y=grpmock.Zfof;
for i=1:numel(massbin)-1
    ff=f&grpmock.Luminosity>10.^massbin(i)&grpmock.Luminosity<=10.^massbin(i+1);
    semilogx(x(ff),y(ff),'o','color',colors(i,:),'markersize',(massbin(i)-5)^1.4);hold on;
end
%%
% grpmock.LuminosityRaw=grpmock.Luminosity./(Bc+Bn./sqrt(grpmock.Mult)+Bz./sqrt(grpmock.Zfof));
figure;
x=grpmock.Zfof;y=grpmock.Luminosity./grpmock.Mult;
f=grpmock.Mult>2;
massbin=[2,4,6,10,20,50,100];
% massbin=9:0.5:13;
colors=colormap(hsv(numel(massbin)));
for i=1:numel(massbin)-1
    ff=f&grpmock.Mult>massbin(i)&grpmock.Mult<=massbin(i+1);
    semilogy(x(ff),y(ff),'o','color',colors(i,:),'markersize',sqrt(massbin(i)));hold on;
end
%%
f=grpmock.Mult>2;
x=grpmock.Zfof(f);y=grpmock.LuminosityRaw./grpmock.Mult;
% cftool(x(f),log10(1./y(f)));
figure;loglog(1./y(f),ffz(x),'.');hold on;
plot(grp.Mult(grp.Mult>2)./grp.LuminosityRaw(grp.Mult>2),ffz(grp.Zfof(grp.Mult>2)),'y.');
plot([1e-14,1e-7],[1e-14,1e-7],'r');
%% fitted (N/L)(z)
h=0.7;
dm=@(z) 5*log10(lum_dist(0.3,z)*1e6/h/10);
pfz=[-103.5,64.12,-16.84,-8.77];
fz=@(z) polyval(pfz,z);
% sfz=@(x) (10*exp(-x))./(3*x.^(0.3)) - (10*gammainc(x,0.7,'upper')*gamma(0.7))/3;
% sfz=@(x) (5*exp(-x))./x.^(1/5) - 5*gammainc(x,4/5,'upper')*gamma(0.8);
sfz=@(x,a) -(exp(-x)./(a+1).*x.^(a+1))+gammainc(x,2+a,'upper')*gamma(2+a)/(a+1);
% sfz=@(l) quadgk(@(x) x.^(-1.2).*exp(-x),l,inf)
ffz=@(z,ml,Mc,a) sfz(10.^(-0.4*(ml-Mc-dm(z))),a)./10.^(-0.4*(Mc-4.76))/0.7^2/gamma(2+a);
ffz2=@(z,ml,Mc,a) sfz(10.^(-0.4*(ml-Mc-dm(z))),a)./sfz(10.^(-0.4*(-13-Mc)),a); %corrected to M=-13, slightly fainter than the limiting m=19.4 at z=0.01 is M=-13.78
%%
ml=[19.4,19.4,19.4]; %groups optimized for ml=19.4, even gama-12 has real ml=19.8, but we used 19.4 version.
Mc=-20.73+5*log10(0.7);
a=-1.26;
grp.IMult=grp.Mult./ffz2(grp.Zfof,19.4,Mc-0.7*grp.Zfof,a);
figure;
loglog(grp.Mult,grp.IMult,'.');hold on;
plot([1,1e3],[1,1e3],'k-');
%%
figure;
colors=colormap(hsv(5));
for sky=1:3
    grpcat{sky}.MultPred=grpcat{sky}.LuminosityRaw.*ffz(grpcat{sky}.Zfof,ml(sky),Mc-0.7*grpcat{sky}.Zfof,a);
    loglog(grpcat{sky}.Mult,grpcat{sky}.MultPred,'.','color',colors(sky,:));hold on;
end
plot([1,100],[1,100],'k')
%%
z=0.01:0.05:0.5;
cftool(z,log(ffz(z,19.4,Mc-0.7*z,a)));
%% fit ffz(z) with exponentional decay
figure;
loglog(grp.Mult,(grp.Luminosity.*exp(-21.5-14.48*grp.Zfof))./grp.Mult,'.');hold on;
plot([1,500],[1,1],'k')
xlabel('N');ylabel('Npred/N');
%%
figure;
z=0.1:0.1:0.5;
plot(z,fz(z)-log10(ffz(z)),'r');
z=0.1:0.1:0.5;
plot(z,fz(z)-log10(ffz(z)),'r');
hold on;
% plot(z,log10(ffz(z)),'g');
%%
figure;
z=0:0.02:1;
semilogy(z,ffz(z),'.')
hold on;
% plot(z,log10(ffz(z)),'g');
%%
figure;
z=0:0.02:1;
semilogy(z,ffz(z),'.')
%% Nabs dependence in the raw mock, queried from the DGalaxies database, Bower2006a.
data=importdata('Mhalo-Ngal-FixLgrp.mock',',',26);
data=data.data;
data(:,1)=data(:,1)*8.6e8;%Msun/h
mbrd=log10(minmax(data(:,3)'));
nbin=logspace(mbrd(1),mbrd(2),40);
colors=spring(numel(nbin));
myfigure;
for i=1:numel(nbin)-1
    f=data(:,3)>=nbin(i)&data(:,3)<nbin(i+1);
loglog(data(f,2),data(f,1),'o','color',colors(i,:));
hold on;
end
xscale('linear');
yscale('log');
xlabel('$L_{grp} [\rm{L}_\odot/h^2]$');
ylabel('$M_{h} [\rm{M}_\odot/h]$');
xlim([0.8e11,2.2e11]);
colormap(spring(numel(nbin)));%caxis([1, numel(nbin)]);
hcb = colorbar('YTick',(log10([10,100,1000])-mbrd(1))/(mbrd(2)-mbrd(1))*(numel(nbin)-1)+1,'YTickLabel',...
{'10','100','1000'});
text(2.35e11, 1e10,'$N_{abs}$','fontsize',18)
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/Mhalo-Lhalo-N.eps');
%%
myfigure;
% mbrd=(minmax(data(:,2)'));
% nbin=linspace(mbrd(1),mbrd(2),30);
% colors=colormap(jet(numel(nbin)));
% for i=1:numel(nbin)-1
%     f=data(:,2)>=nbin(i)&data(:,2)<nbin(i+1);
% loglog(data(f,3),data(f,1),'o','color',colors(i,:));
% hold on;
% end
loglog(data(:,3),data(:,1),'s');
xscale('log');
yscale('log');
xlabel('$N_{gal}$');
ylabel('$M_{halo} [\rm{M}_\odot/h]$');
% xlim([0.85e11,1.05e11]);
% colormap(jet(numel(nbin)));%caxis([1, numel(nbin)]);
% hcb = colorbar('YTick',(log10([1,10,100,1000])-mbrd(1))/(mbrd(2)-mbrd(1))*(numel(nbin)-1)+1,'YTickLabel',...
% {'N','10','100','1000'});
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3/Mhalo-Ngal-FixLhalo.eps');
%% M=A(L)
data=[3.0851e+10   1.08705e+10  12.22      0.5089                  
7.16576e+10   1.39753e+10  13.01      0.2178                  
1.57631e+11   4.07415e+10  13.45      0.1417                  
3.79365e+11   9.89501e+10  13.71      0.2087                  
9.70493e+11   3.28707e+11  14.21      0.2076                  
2.90161e+12   6.07367e+11  14.79      0.3664];
xbrd=[1e10,5e10; 5e10,1e11; 1e11,2.5e11 ;2.5e11,6e11 ;6e11,2e12 ;2e12,5e12];
myfigure;
for vol=1
    load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
f=grpmock.Mult>2;
x=logspace(9,13,12);
[xmed,ymed,yl]=skeleton(grpmock.Luminosity(f),(grpmock.MIter(f)),x,0.68);
h3=plot(xmed,ymed,'g');hold on;
end
h33=plot(xmed,yl(:,1),'g--');
plot(xmed,yl(:,2),'g--');
x=data(:,1);
data(:,3)=data(:,3)-0.07; %correct for sysmatics
h1=ploterr(x,10.^data(:,3),{xbrd(:,1),xbrd(:,2)},{10.^(data(:,3)-data(:,4)),10.^(data(:,3)+data(:,4))},'ro','logxy');hold on;
x=logspace(10,13);
% h2=plot(x,10.^(14.23+1.09*log10(x/1e12)),'k');
h22=plot(x,10.^(14.23-0.07+1.1*log10(x/1e12)),'k-');
% h3=plot(x,14.18+1.051*log10(x/1e12),'k--'); %without Xoffset cut
xlabel('Group Luminosity[L$_\odot/h^2$]');
ylabel('Halo Mass[M$_\odot/h$]');
xlim([8e9,6e12]);ylim([2e11,2e15]);
set(gca,'xscale','log','yscale','log');
l=legend([h1(1),h22,h3,h33],'MLWL-Binned','MLWL-Parametric','Mock-Median','Mock-$1\sigma$ Bounds');set(l,'location','southeast','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/extra/M-L-poster.eps');