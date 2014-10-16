clear
clc
load G3Cv4up8/mockcat_1.mat
load G3Cv4up8/G3Cv4up8.mat
%% LumMass 
data=[
   7.44696e+12   1.41114e+12  12.13      0.7874                  
1.77823e+13   5.56486e+12  13.06      0.1747                  
4.55094e+13   1.13759e+13  13.45      0.1505
1.12747e+14   3.32744e+13  13.79      0.1747                  
3.04851e+14   9.30861e+13  14.02      0.3086
]; 
xbrd=[5e12,1e13; 1e13,3e13; 3e13,7e13; 7e13,2e14; 2e14,6e14];
myfigure;
for vol=1:9
    load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
f=grpmock.Mult>2;
x=logspace(11,15,12);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.LumMass(f),log10(grpmock.MIter(f)),x,0.68);
h3=plot(xmed,ymed,'g');hold on;
end
plot(xmed,yl,'g--');
x=data(:,1);
h1=ploterr(x,data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro','logx');hold on;
x=logspace(10,13);
% h2=plot(x,14.24+1.09*log10(x/1e12),'k');
xlabel('LumMass');
ylabel('$\log(M_h$[M$_\odot/h$])');
% xlim([8e9,6e12]);
set(gca,'xscale','log');
l=legend([h1(1),h2,h3],'Binned','Global','Mock');set(l,'location','southeast','interpreter','latex');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/M-L.eps');
%% LumBias, linear scale fit
data=[7.44696e+12   1.41114e+12  0.2076     0.3433                  
1.77823e+13   5.56486e+12  0.6283     0.2657                  
4.55094e+13   1.13759e+13  0.6282     0.2193                  
1.12747e+14   3.32744e+13  0.6091     0.233                   
3.04851e+14   9.30861e+13  0.3654     0.261
9.37345e+14   2.52402e+14  0.3065     0.3605
];
xbrd=[5e12,1e13; 1e13,3e13; 3e13,7e13; 7e13,2e14; 2e14,6e14; 6e14,2e15];
myfigure;
for vol=1:9
    load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
f=grpmock.Mult>2;
x=logspace(11.8,16,10);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.LumMass(f),(grpmock.MIter(f)./grpmock.LumMass(f)),x,0.68);
h3=plot(xmed,ymed,'g');hold on;
end
plot(xmed,yl,'g--');
x=data(:,1);
h1=ploterr(x,data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro','logx');hold on;
x=logspace(12,16,2);
% h2=plot(x,14.24+1.09*log10(x/1e12),'k');
xlabel('$L_{grp}$[L$_\odot/h^2$]');
ylabel('$\log(M_h$[M$_\odot/h$])');
xlim([1e12,1e16]);
ylim([-1.5,1]);
set(gca,'xscale','log');
l=legend([h1(1),h2,h22,h3],'Binned','Global','Systematic','Mock');set(l,'location','southeast','interpreter','latex');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/M-L.eps');
%% LumBias, log scale fit
data=[7.44696e+12   1.41114e+12 -0.6815     0.7162                  
1.77823e+13   5.56486e+12 -0.2017     0.1835                  
4.55094e+13   1.13759e+13 -0.2017     0.1516                  
1.12747e+14   3.32744e+13 -0.2153     0.1661                  
3.04851e+14   9.30861e+13 -0.4368     0.31                    
9.37345e+14   2.52402e+14 -0.5122     0.5094 
];
xbrd=[5e12,1e13; 1e13,3e13; 3e13,7e13; 7e13,2e14; 2e14,6e14; 6e14,2e15];
myfigure;
for vol=1:9
    load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
f=grpmock.Mult>2;
x=logspace(11.8,16,10);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.LumMass(f),log10(grpmock.MIter(f)./grpmock.LumMass(f)),x,0.68);
h3=plot(xmed,ymed,'g');hold on;
end
plot(xmed,yl,'g--');
x=data(:,1);
h1=ploterr(x,data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro','logx');hold on;
x=xbrd([1,end]);
h2=plot(x,[-0.28,-0.28],'k-');
h22=plot(x,[-0.28,-0.28]-0.09,'k--');
% h2=plot(x,-0.12+0.01*log10(x/1e12),'k:'); %N>=10.
xlabel('$M_{lum}$[M$_\odot/h$]');
ylabel('$\log(M_h/M_{lum})$');
% xlim([1e12,1e16]);
xlim([4e11,3e15]);
ylim([-1.5,1]);
set(gca,'xscale','log');
l=legend([h1(1),h2,h22,h3],'Binned','Global','Systematic','Mock');set(l,'location','northeast','interpreter','latex');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/LumBias.eps');
%% DynBias, log scale fit
data=[4.51771e+12   2.74885e+12 -0.6248     0.792                   
1.86774e+13   5.8161e+12 -0.6095     0.3526                  
4.66888e+13   1.12103e+13 -0.2908     0.1651                  
1.18891e+14   3.58124e+13 -0.3678     0.138                   
3.2316e+14   1.02288e+14 -1.58       0.9296                  
9.34711e+14   2.96602e+14 -0.5167     0.2546
];
xbrd=[5e11,1e13; 1e13,3e13; 3e13,7e13; 7e13,2e14; 2e14,6e14; 6e14,2e15];
myfigure;
for vol=1:9
    load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
f=grpmock.Mult>2;
x=logspace(11,16,10);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.DynMass(f),log10(grpmock.MIter(f)./grpmock.DynMass(f)),x,0.68);
h3=plot(xmed,ymed,'g');hold on;
end
plot(xmed,yl,'g--');
x=data(:,1);
h1=ploterr(x,data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro','logx');hold on;
x=xbrd([1,end]);
h2=plot(x,0.08-0.31*log10(x/1e12),'k-');
h22=plot(x,0.08-0.3-(0.31+0.04)*log10(x/1e12),'k--');
% h2=plot(x,14.24+1.09*log10(x/1e12),'k');
xlabel('$M_{dyn}$[M$_\odot/h$]');
ylabel('$\log(M_h/M_{dyn})$');
xlim([4e11,3e15]);
ylim([-4,2]);
set(gca,'xscale','log');
set(gca,'xtick',[1e12,1e13,1e14,1e15]);
l=legend([h1(1),h2,h22,h3],'Binned','Global','Systematic','Mock');set(l,'location','northeast','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/DynBias.eps');
%% LumBias, mult dependence
data=[3.282903   0.450410 -0.3616     0.1909                  
5.622797   0.772559 -0.5123     0.2444                  
9.495979   1.347467 -0.2595     0.1838                  
15.000000   1.681543  0.07959    0.1988                  
21.324324   1.973528  0.07684    0.2559                  
31.457143   4.842330 -0.1882     0.3459                                 
];
data=[4.936483   5.587022 -0.2792     0.09162                 
8.632300   8.986510 -0.2586     0.1047                  
11.979710   11.824423 -0.148      0.1068                  
14.855792   14.371482 -0.136      0.1181                  
18.344961   17.528968  0.0008766    0.1122                  
22.368098   21.029619 -0.05237    0.1396                  
26.294643   24.373814 -0.08884    0.1639                  
33.816667   31.389219 -0.04862    0.178]; %cumulative
xmin=[2 4 6 8 10 12 15 20];
% xbrd=[2,4; 4,7; 7,12; 12,18; 18,25; 25,50];
myfigure;
for vol=1:9
    load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
f=grpmock.Mult>2;
x=[xbrd(:,1);xbrd(end)];
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.Mult(f),log10(grpmock.MIter(f)./grpmock.LumMass(f)),x+1,0.68);
h3=plot(xmed,ymed,'g');hold on;
end
plot(xmed,yl,'g--');
x=data(:,1);
h1=ploterr(x,data(:,3),{xbrd(:,1)+1,xbrd(:,2)},data(:,4),'ro','logx');hold on;
x=xbrd([1,end]);
h2=plot(x,[-0.28,-0.28],'k-');
h22=plot(x,[-0.28,-0.28]-0.09,'k--');
% h2=plot(x,-0.12+0.01*log10(x/1e12),'k:'); %N>=10.
xlabel('$M_{lum}$[L$_\odot/h^2$]');
ylabel('$\log(M_h/M_{lum})$');
% xlim([1e12,1e16]);
ylim([-1.5,1]);
set(gca,'xscale','log');
l=legend([h1(1),h2,h22,h3],'Binned','Global','Systematic','Mock');set(l,'location','northeast','interpreter','latex');
%% LumBias, mult dependence
data=[4.936483   5.587022 -0.2792     0.09162                 
6.837898   7.386943 -0.224      0.09282                 
8.632300   8.986510 -0.2586     0.1047                  
10.524123   10.600370 -0.1434     0.09945                 
11.979710   11.824423 -0.148      0.1068                  
13.410448   13.069683 -0.1428     0.1121                  
14.855792   14.371482 -0.136      0.1181                  
16.645062   15.999053 -0.09986    0.1217                  
18.344961   17.528968  0.0008766    0.1122                  
20.243902   19.213284  0.004493    0.119                   
22.368098   21.029619 -0.05237    0.1396                  
23.907143   22.318441 -0.05932    0.1484                  
25.368852   23.558160 -0.0678     0.1547                  
26.294643   24.373814 -0.08884    0.1639                  
30.802632   28.495161 -0.1324     0.1937                  
33.816667   31.389219 -0.04862    0.178]; %cumulative
datadyn=[4.936483   5.587022 -0.6019     0.1045                  
6.837898   7.386943 -0.4301     0.09731                 
8.632300   8.986510 -0.3803     0.1021                  
10.524123   10.600370 -0.2526     0.09668                 
11.979710   11.824423 -0.2137     0.09973                 
13.410448   13.069683 -0.2025     0.1047                  
14.855792   14.371482 -0.2093     0.1128                  
16.645062   15.999053 -0.1832     0.1175                  
18.344961   17.528968 -0.06323    0.1063                  
20.243902   19.213284 -0.08247    0.117                   
22.368098   21.029619 -0.1229     0.1354                  
23.907143   22.318441 -0.1151     0.1411                  
25.368852   23.558160 -0.1366     0.1496                  
26.294643   24.373814 -0.1467     0.155                   
30.802632   28.495161 -0.1756     0.1798                  
33.816667   31.389219 -0.1335     0.1757]; %cumulative
xmin=[2:15,18,20];
% xbrd=[2,4; 4,7; 7,12; 12,18; 18,25; 25,50];
myfigure;
for vol=1
    load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
    y=log10(grpmock.MIter./grpmock.LumMass);
    ymed=xmin;
    yl=[xmin;xmin];
    for i=1:numel(xmin)
        ymed(i)=median(y(grpmock.Mult>xmin(i)));
        yl(:,i)=prctile(y(grpmock.Mult>xmin(i)),[0.5-0.683/2,0.5+0.683/2]*100);
    end
% h3=plot(xmin+1,10.^ymed,'g');hold on;
end
% plot(xmin+1,10.^yl,'g--');
h1=ploterr(xmin+1,10.^data(:,3),[],{10.^(data(:,3)+data(:,4)),10.^(data(:,3)-data(:,4))},'ro');hold on;
% h1=plot(xmin+1,10.^data(:,3),'r-');hold on;
% h11=plot(xmin+1,[10.^(data(:,3)+data(:,4)),10.^(data(:,3)-data(:,4))],'r--');hold on;
h2=plot(xmin+1,10.^datadyn(:,3),'b-');hold on;
h22=plot(xmin+1,[10.^(datadyn(:,3)+datadyn(:,4)),10.^(datadyn(:,3)-datadyn(:,4))],'b--');hold on;

plot([2,22],[1,1],'k--');
xlabel('$N_{min}$');
ylabel('$M_h/M_{G3C}$');
% xlim([1e12,1e16]);
% ylim([-1.5,1]);
% set(gca,'xscale','log');
l=legend([h1(1),h2],'LumMass','DynMass');set(l,'location','southeast','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/G3CBiasMult.eps');
%% DynBias, mult dependence
data=[4.936483   5.587022 -0.6019     0.1045                  
6.837898   7.386943 -0.4301     0.09731                 
8.632300   8.986510 -0.3803     0.1021                  
10.524123   10.600370 -0.2526     0.09668                 
11.979710   11.824423 -0.2137     0.09973                 
13.410448   13.069683 -0.2025     0.1047                  
14.855792   14.371482 -0.2093     0.1128                  
16.645062   15.999053 -0.1832     0.1175                  
18.344961   17.528968 -0.06323    0.1063                  
20.243902   19.213284 -0.08247    0.117                   
22.368098   21.029619 -0.1229     0.1354                  
23.907143   22.318441 -0.1151     0.1411                  
25.368852   23.558160 -0.1366     0.1496                  
26.294643   24.373814 -0.1467     0.155                   
30.802632   28.495161 -0.1756     0.1798                  
33.816667   31.389219 -0.1335     0.1757]; %cumulative
datalum=[4.936483   5.587022 -0.2792     0.09162                 
6.837898   7.386943 -0.224      0.09282                 
8.632300   8.986510 -0.2586     0.1047                  
10.524123   10.600370 -0.1434     0.09945                 
11.979710   11.824423 -0.148      0.1068                  
13.410448   13.069683 -0.1428     0.1121                  
14.855792   14.371482 -0.136      0.1181                  
16.645062   15.999053 -0.09986    0.1217                  
18.344961   17.528968  0.0008766    0.1122                  
20.243902   19.213284  0.004493    0.119                   
22.368098   21.029619 -0.05237    0.1396                  
23.907143   22.318441 -0.05932    0.1484                  
25.368852   23.558160 -0.0678     0.1547                  
26.294643   24.373814 -0.08884    0.1639                  
30.802632   28.495161 -0.1324     0.1937                  
33.816667   31.389219 -0.04862    0.178]; %cumulative
xmin=[2:15,18,20];
% xbrd=[2,4; 4,7; 7,12; 12,18; 18,25; 25,50];
myfigure;
for vol=1
    load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
    y=log10(grpmock.MIter./grpmock.DynMass);
    ymed=xmin;
    yl=[xmin;xmin];
    for i=1:numel(xmin)
        ymed(i)=median(y(grpmock.Mult>xmin(i)));
        yl(:,i)=prctile(y(grpmock.Mult>xmin(i)),[0.5-0.683/2,0.5+0.683/2]*100);
    end
h3=plot(xmin+1,10.^ymed,'g');hold on;
end
plot(xmin+1,10.^yl,'g--');
h1=ploterr(xmin+1,10.^data(:,3),[],{10.^(data(:,3)+data(:,4)),10.^(data(:,3)-data(:,4))},'ro-','logx');hold on;
% h1=plot(xmin+1,10.^data(:,3),'r-');hold on;
% h1=plot(xmin+1,[10.^(data(:,3)+data(:,4)),10.^(data(:,3)-data(:,4))],'r--');hold on;
h1=plot(xmin+1,10.^datalum(:,3),'b-');hold on;
plot([2,22],[1,1],'k--');
xlabel('$N_{min}$');
ylabel('$M_h/M_{dyn}$');
% xlim([1e12,1e16]);
ylim([0,2]);
% set(gca,'xscale','log');
l=legend([h1(1),h3],'Binned','Mock');set(l,'location','northeast','interpreter','latex');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/DynBiasMult.eps');