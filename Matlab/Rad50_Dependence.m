clear
clc
load G3Cv4up8/mockcat_1.mat
load G3Cv4up8/G3Cv4up8.mat
%% M=A(R)
xbrd=[0.02,0.08; 0.08,0.24; 0.24,0.5; 0.5,2];
data=[
0.0239661   0.00777136  11.55      2.139                   
0.0547403   0.0101621  12.48      0.582                   
0.107018   0.0196825  12.99      0.2348                  
0.19745   0.0360921  13.25      0.1529                  
0.357995   0.067244  13.13      0.3163                  
0.649148   0.124558  13.68      0.4731];
xbrd=[0.01,0.037; 0.037,0.072; 0.072,0.14; 0.14,0.27; 0.27,0.52; 0.52,1.1];
% data=[0.0145983   0.00290917  11.66      2.664                   
% 0.0282677   0.0050438  11.61      2.793                   
% 0.0547403   0.0101621  12.58      0.5063                  
% 0.107018   0.0196825  12.96      0.2436                  
% 0.19745   0.0360921  13.24      0.1521                  
% 0.357995   0.067244  13.03      0.3384                  
% 0.649148   0.124558  13.45      0.6817]; %without Xoffset cut
% xbrd=[0.01,0.02; 0.02,0.037; 0.037,0.072; 0.072,0.14; 0.14,0.27;
% 0.27,0.52; 0.52,1.1];
myfigure;
for vol=1:9
     load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
f=grpmock.Mult>2;
x=grpmock.Rad50(f);y=log10(grpmock.MIter(f));
xbin=logspace(-3,0.3,10);
[xmed,ymed,yl]=skeleton(x,y,xbin,0.68);
h2=plot(xmed,ymed,'g');hold on;
end
plot(xmed,yl,'g--');
h1=ploterr(data(:,1),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro','logx');hold on;
x=[0.01,1];
h3=plot(x,13.85+0.98*log10(x),'k-');
h33=plot(x,13.85-0.03+(0.98+0.05)*log10(x),'k--');
xlabel('$R_{50}$[(Mpc$/h)$]');
ylabel('$\log(M_h[\rm{M}_{\odot}/h])$');
set(gca,'xscale','log');
plot([5e-3,5e-3],[9,15],'k--');
%  xlim([90,2100]);
% ylim([11,15.5]);
l=legend([h1(1),h3,h33,h2],'Binned','Global','Systematic','Mock');set(l,'location','southeast','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/M-R.eps');
%% M=A(AH)
data=[0.133783   0.0484982  12.8       1.03                    
0.928755   0.483627  12.65      0.2826                  
4.45234   2.10506  13.31      0.1604                  
18.1153   6.8864  13.99      0.1409                  
62.1664   17.0101  14.67      0.1765];
xbrd=[0.005,0.2; 0.2,2; 2,10; 10,40; 40,222];
myfigure;
% for i=1:9
%     grpmock=mock{i}.grpmock;
f=grpmock.Mult>2;
x=grpmock.d3area(f);y=log10(grpmock.MIter(f));
xbin=logspace(-4,3,8);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(x,y,xbin,0.68);
h2=plot(xmed,ymed,'g');hold on;
plot(xmed,yl(:,1),'g--');
plot(xmed,yl(:,2),'g--');
% end
h1=ploterr(data(:,1),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro','logx');hold on;
x=data(:,1);
% h3=plot(x,13.87+1.23*log10(x),'k-');
xlabel('$V_{N}$[(Mpc$/h)^3$]');
ylabel('$\log(M_h[\rm{M}_{\odot}/h])$');
set(gca,'xscale','log');
%  xlim([90,2100]);
% ylim([11,15.5]);
l=legend([h1(1),h3,h2],'Binned','Global','Mock');set(l,'location','southeast','interpreter','latex');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/v3/M-Vol.eps');
%% M=A(z)
data=[0.0678874   0.0206498  11.9       1.401                   
0.150748   0.0281476  13.09      0.1497                  
0.248102   0.0300998  13.29      0.2135                  
0.332109   0.0256205  13.25      0.6272                  
0.414707   0.0149001  7.026      nan];
x=data(:,1);
xbrd=[0,0.1;0.1,0.2;0.2,0.3;0.3,0.4;0.4,0.5];
myfigure;
ploterr(x,data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'o');
hold on;
plot([0,0.5],[8.92,8.92],'r-')
xlabel('Redshift');
ylabel('$\log(M)$');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/M-z.eps');
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
%% M=A(R)
data=[
0.0239661   0.00777136  11.55      2.139                   
0.0547403   0.0101621  12.48      0.582                   
0.107018   0.0196825  12.99      0.2348                  
0.19745   0.0360921  13.25      0.1529                  
0.357995   0.067244  13.13      0.3163                  
0.649148   0.124558  13.68      0.4731];
xbrd=[0.01,0.037; 0.037,0.072; 0.072,0.14; 0.14,0.27; 0.27,0.52; 0.52,1.1];
myfigure;
for vol=1
     load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
f=grpmock.Mult>2;
x=grpmock.Rad50(f);y=(grpmock.MIter(f));
xbin=logspace(-3,0.3,10);
[xmed,ymed,yl]=skeleton(x,y,xbin,0.68);
h2=plot(xmed,ymed,'g');hold on;
end
h22=plot(xmed,yl,'g--');
data(:,3)=data(:,3)-0.03;
h1=ploterr(data(:,1),10.^data(:,3),{xbrd(:,1),xbrd(:,2)},{10.^(data(:,3)-data(:,4)),10.^(data(:,3)+data(:,4))},'ro','logx');hold on;
x=[0.01,1];
h33=plot(x,10.^(13.85-0.03+(0.98+0.05)*log10(x)),'k-');
xlabel('Group Radius[(Mpc$/h)$]');
ylabel('Halo Mass$[\rm{M}_{\odot}/h]$');
set(gca,'xscale','log','yscale','log');
% plot([5e-3,5e-3],10.^[9,15],'k--');
xlim([5e-3,10]);
% ylim([11,15.5]);
l=legend([h1(1),h33,h2,h22(1)],'MLWL-Binned','MLWL-Parametric','Mock-Median','Mock-$1\sigma$ Bounds');set(l,'location','southeast','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/extra/M-R-poster.eps');