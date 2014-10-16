clear
clc
load G3Cv4up8/mockcat_1.mat
load G3Cv4up8/G3Cv4up8.mat
%%
% newgal=fits_load_bintable('201308/TilingCatv40_kcorr_z00v03.fits');
tmp=load('201310/mockcat_1.mat','galmock'); newgalmock=tmp.galmock; clear tmp;
%% adopting GAMA-II dN/dV
Nav=load('GAMAII-Ndens-z.txt','-ascii');
avNdens=@(z) pchip(Nav(:,1),Nav(:,2),z);
grpmock.VolMult2=grpmock.Mult./avNdens(grpmock.Zfof);
%% M=A(V)
data=[
181.758   59.369  12.91      0.1816                  
575.397   186.928  13.29      0.1519                  
1798.16   649.564  13.78      0.2                     
5769.01   2058.02  14.13      0.4463                  
13834.1   1347.12  15.19      0.6859]; %N>2 
xbrd=[90,310; 310,1050; 1050,3600; 3600,12000; 12000,42000];
myfigure;
for vol=1:9
    load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
f=grpmock.Mult>2;
x=grpmock.VolMult(f);y=log10(grpmock.MIter(f));
[xmed,ymed,yl,xm,ym,yerr]=skeleton(x,y,logspace(1.3,5,15),0.68);
h2=plot(xmed,ymed,'g');hold on;
end
plot(xmed,yl(:,1),'g--');
plot(xmed,yl(:,2),'g--');
h1=ploterr(data(:,1),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro','logx');hold on;
h3=plot(xmed,10.49+1.03*log10(xmed),'k-');
h33=plot(xmed,10.49-0.1+1.04*log10(xmed),'k--');
xlabel('$V_{N}$[(Mpc$/h)^3$]');
ylabel('$\log(M_h[\rm{M}_{\odot}/h])$');
set(gca,'xscale','log');
%  xlim([90,2100]);
% ylim([11,15.5]);
l=legend([h1(1),h3,h33,h2],'Binned','Global','Systematic','Mock');set(l,'location','southeast','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/M-Vol.eps');
%% M=A(1+z)
data=[
    %0.0470461   0.0137519  5.367      453.7                   
0.109357   0.021428  12.65      0.2882                  
0.17942   0.0219766  13.34      0.1377                  
0.265169   0.0242822  13.11      0.3475                  
0.339685   0.0239398  13.64      0.4202                  
0.414394   0.0148172  13.48      2.063 ];
x=data(:,1);
xbrd=1+[0.07,0.14; 0.14,0.22; 0.22,0.31; 0.31,0.40; 0.40,0.5];
myfigure;
for vol=1:9
    load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
f=grpmock.Mult>2;
x=grpmock.Zfof(f)+1;y=log10(grpmock.MIter(f));
[xmed,ymed,yl,xm,ym,yerr]=skeleton(x,y,logspace(0,log10(1.5),15),0.68);
h2=plot(xmed,ymed,'g');hold on;
end
plot(xmed,yl,'g--');
h1=ploterr(1+data(:,1),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro');
x=[1.,1.5];
h3=plot(x,12.38+9.7*log10(x),'k-');
h3=plot(x,12.46+8.395*log10(x),'r--'); %after excluding z>1.0 sources.
xlabel('1+z');
ylabel('$\log(M_h[\rm{M}_{\odot}/h])$');
l=legend([h1(1),h3,h2],'Binned','Global','Mock');set(l,'location','northwest','interpreter','latex');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/M-z.eps');
%% M=A(z)
data=[
    %0.0470461   0.0137519  5.367      453.7                   
0.109357   0.021428  12.65      0.2882                  
0.17942   0.0219766  13.34      0.1377                  
0.265169   0.0242822  13.11      0.3475                  
0.339685   0.0239398  13.64      0.4202                  
0.414394   0.0148172  13.48      2.063 ];
xbrd=[0.07,0.14; 0.14,0.22; 0.22,0.31; 0.31,0.40; 0.40,0.5];
myfigure;
for vol=1:9
    load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
f=grpmock.Mult>2;
x=grpmock.Zfof(f);y=log10(grpmock.MIter(f));
[xmed,ymed,yl,xm,ym,yerr]=skeleton(x,y,logspace(-2,log10(0.52),15),0.68);
h2=plot(xmed,ymed,'g');hold on;
end
plot(xmed,yl,'g--');
ploterr(data(:,1),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'o');
x=logspace(-1.2,log10(0.5),10);
plot(x,14.24+1.46*log10(x),'k-');
xlabel('z');
ylabel('$\log(M)$');
xscale('log');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/M-z.eps');
%% M-sigma
data=[141.452   22.5981  8.988      0.3632                  
249.279   43.1649  8.766      0.2191                  
429.424   72.3064  8.928      0.1021                  
746.285   119.504  9.066      0.1368                  
1295.35   166.13  7.671      nan ]; %M=A*f(N)*sigma^2
xbrd=[0.1e3,0.18e3; 0.18e3,0.33e3; 0.33e3,0.60e3; 0.60e3,1.10e3; 1.10e3,2.0e3];
myfigure;
h=ploterr(data(:,1),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'o','logx');hold on;
plot([100,2000],[8.9,8.9],'k-');
xlabel('$\sigma_v$[km/s]');
ylabel('A');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/M-sigma2.eps');
%% M=A(sigma)
data=[141.452   22.5981  11.66      1.583                   
249.279   43.1649  12.93      0.1919                  
429.424   72.3064  13.33      0.1744                  
746.285   119.504  14.24      0.1548                  
1295.35   166.13  14.14      1.361];
xbrd=[0.1e3,0.18e3; 0.18e3,0.33e3; 0.33e3,0.60e3; 0.60e3,1.10e3; 1.10e3,2.0e3];
myfigure;
mfit=10.^8.9*min(grp.Mult/14,1.0).^1.8.*grp.VelDisp.^2; %MultFit
semilogx(grp.VelDisp(grp.Mult>2),log10(mfit(grp.Mult>2)),'o','color',[0.85,0.85,0.85]);hold on;
semilogx(grp.VelDisp(grp.Mult==2),log10(mfit(grp.Mult==2)),'s','color',[0.95,0.85,0.85]);hold on;
h=ploterr(data(:,1),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'o','logx');hold on;
% plot([100,2000],[8.9,8.9],'k-');
f=grpmock.Mult>2;
multbin=[2,4,6,10,20,50,100];
colors=colormap(cool(numel(multbin)));
x=grpmock.VelDisp(f);y=log10(grpmock.MIter(f));
for i=1:numel(multbin)-1
    ff=f&grpmock.Mult>multbin(i)&grpmock.Mult<=multbin(i+1);
% plot(x(ff),y(ff),'.','color',colors(i,:));hold on;
end
% plot(grpmock.VelDisp(f),log10(mockmass.MIter(f)),'.')
x=logspace(log10(80),log10(2000),10);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.VelDisp(f),log10(mockmass.MIter(f)),x,0.68);
plot(xmed,ymed,'r');hold on;
plot(xmed,yl(:,1),'r--');
plot(xmed,yl(:,2),'r--');
xlabel('$\sigma_v$[km/s]');
ylabel('M');
set(gca,'xscale','log');
% xlim([80,2100]);
% ylim([11,15.5]);
% print('-depsc','/work/Projects/Lensing/outputv4/paper/M-sigma.eps');
%%
data=[3.282903   0.450410  7.802      0.2497                  
5.622797   0.772559  8.206      0.2295                  
9.495979   1.347467  8.658      0.1398                  
15.000000   1.681543  8.906      0.181                   
21.324324   1.973528  9.072      0.2084                  
29.933333   2.874408  8.735      0.2653                  
48.750000   11.255554  8.93       0.3262 ];

myfigure;
x=data(:,1);
xbrd=[3,4; 5,7; 8,12; 13,18; 19,25; 26,35; 36,100];
% x=log10(x);
% xbrd=log10(xbrd);
h=ploterr(x,data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'o','logx');hold on;
x=3:0.1:100;
y=8.92+1.794*log10(min(x/13.64,1));
h1=plot(x,y,'k-');
y=7.8+0.7*log10(x);
h2=plot(x,y,'k--');
f=grpmock.Mult>2;
c=log10(grpmock.VelDisp(f).^2);
c=round(1+(20-1)*(c-min(c))/(max(c)-min(c)));
colors=colormap(cool(20));
x=grpmock.Mult(f);y=log10(mockmass.MIter(f)./grpmock.VelDisp(f).^2);
for i=1:numel(x)
plot(x(i),y(i),'.','color',colors(c(i),:));hold on;
end
% plot(grpmock.Mult(f),log10(mockmass.MIter(f)./grpmock.VelDisp(f).^2),'y.')
% x=[xbrd(:,1)',xbrd(end)];
x=[2.5:1.5:12,15,20,40,70,100];
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.Mult(f),log10(mockmass.MIter(f)./grpmock.VelDisp(f).^2),x,0.68);
h3=plot(xmed,ymed,'r');hold on;
plot(xmed,yl(:,1),'r--');
plot(xmed,yl(:,2),'r--');
xlabel('Multiplicity');
l=legend([h(1),h1,h2,h3],'Data','$A\cdot$min$(N/N_0,1)^b$','$A\cdot N^b$','Mock');
set(l,'interpreter','latex','location','southeast')
ylabel('$\log(M/\sigma_v^2)$');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/M-sigma-N.eps');
%%
myfigure;
mfit=10.^8.9*min(grpmock.Mult/14,1.0).^1.8.*grpmock.VelDisp.^2; %MultFit
% mfit=10.^7.8*grpmock.Mult.^0.7.*grpmock.VelDisp.^2; %MultFit
% mfit=3.7*grpmock.VelDisp.^2.*grpmock.Rad50/4.301e-9; %MultFit
% mfit=10.^2.5*min(grpmock.Mult/14,1.0).*grpmock.TotFluxProxy; %LumFit
% mfit=grpmock.MassProxy;
% mfit=grpmock.LumMass;
f=grpmock.Mult>2&grpmock.Rad50>1e-4;
x=log10(mfit(f));
% y=log10(grpmock.MIter(f));
y=log10(grpmock.MIter(f));
c=log10(min(grpmock.Mult(f)/14,1.0).^1.8);
c=round(1+(20-1)*(c-min(c))/(max(c)-min(c)));
colors=colormap(cool(20));
for i=1:numel(x)
plot(x(i),y(i),'.','color',colors(c(i),:));hold on;
end
[xmed,ymed,yl,xm,ym,yerr]=skeleton(x,y,linspace(10,15,10),0.68);
plot(xmed,ymed,'r');hold on;
plot(xmed,yl(:,1),'r--');
plot(xmed,yl(:,2),'r--');
% ploterr(xm,ym,[],yerr,'bo');
plot([10,15],[10,15],'k-');
xlim([10,16]);ylim([10,16]);
xlabel('$\log(LumMass)$');ylabel('$\log(Mtrue)$');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/extra/MockMass-LumMass.eps');
% semilogx(grpmock.VelDisp(grpmock.Mult==2),log10(mfit(grpmock.Mult==2)),'s','color',[0.95,0.85,0.85]);hold on;
%%
cftool(log10(data(:,1)),data(:,3),[],1./(data(:,4).^2))
% cftool(data(:,1),data(:,3),[])
%% further redshift dependence
data=[0.0678881   0.0205995  8.723      0.2581                  
0.150836   0.0281618  8.965      0.0932                  
0.248039   0.0300626  8.897      0.1749                  
0.332114   0.0256095  9.153      0.2604                  
0.414394   0.0148172  8.436      2.148 ];
x=data(:,1);
xbrd=[0,0.1;0.1,0.2;0.2,0.3;0.3,0.4;0.4,0.5];
myfigure;
ploterr(x,data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'o');
hold on;
plot([0,0.5],[8.92,8.92],'r-')
xlabel('Redshift');
ylabel('$\log(A)$');
print('-depsc','/work/Projects/Lensing/outputv4/paper/M-sigma-z.eps');
%% further Rad50 dependence, comoving
data=[0.052245   0.0269949  9.096      0.4524                   
0.146843   0.028578  8.878      0.1955                  
0.296175   0.076471  8.905      0.121                   
0.660302   0.187962  9.004      0.2909];
x=data(:,1);
xbrd=[0.0,0.1; 0.1,0.2; 0.2,0.5; 0.5,2];
myfigure;
ploterr(x,data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'o');
hold on;
plot([0.,1],[8.92,8.92],'r-')
xlabel('log($R_{cmv}$/[Mpc/h])');
ylabel('$\log(A)$');
% set(gca,'xscale','log');
set(gca,'xlim',[0,1]);
print('-depsc','/work/Projects/Lensing/outputv4/paper/M-sigma-Rad50cmv.eps');
%% further Rad50 dependence, physical
data=[0.0524271   0.0275935  9.121      0.4122                  
0.146155   0.0289109  9.036      0.2748                  
0.285309   0.0725801  8.909      0.1309                  
0.655059   0.152147  8.85       0.5552];
x=data(:,1);
xbrd=[0.0,0.1; 0.1,0.2; 0.2,0.5; 0.5,2];
myfigure;
ploterr(x,data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'o');
hold on;
plot([0.,1],[8.92,8.92],'r-')
xlabel('log($R$/[Mpc/h])');
ylabel('$\log(A)$');
% set(gca,'xscale','log');
set(gca,'xlim',[0,1]);
print('-depsc','/work/Projects/Lensing/outputv4/paper/M-sigma-Rad50.eps');
%%
cftool(log10(data(:,1)),data(:,3),[],1./(data(:,4).^2))
