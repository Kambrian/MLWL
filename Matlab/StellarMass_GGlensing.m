clear
clc
addpath(genpath('/work/Projects/Lensing/code/v8.0/Matlab'))
cd /work/Projects/Lensing/data
load G3Cv4up8/mockcat_1.mat
load G3Cv4up8/G3Cv4up8.mat
galmock.color=galmock.Umodel_abs-galmock.Rmodel_abs;
gal.color=gal.Umodel_abs-gal.Rmodel_abs;
%%
load Eagle
%% raw mocks; dont know what they are.
% tmp=importdata('201310/Bower06.GAMA_20.8.0',' ',6);
% data=tmp.data;
% f=data(:,5)>0&data(:,17)<=19.4;
% rawmock.gid=data(f,3);
% rawmock.rid=data(f,4);
% rawmock.mh=data(f,7);
% rawmock.sm=data(f,30);
% rawmock.sm_disk=data(f,29);
% rawmock.SFR=data(f,26);
% rawmock.z=data(f,10);
%%
load('additional/database/smhm_guo13.mat');
load('additional/database/smhm_bower.mat');
load('additional/database/smhm_font.milli.mat');
guo.Umag=guo.data(:,6);
guo.color=guo.Umag-guo.Rmag;
%%
prepare_smhm_HODs;
%% compare with grp measurement
data=[%7.3685e+09   1.45222e+09  11.51      0.5988                  
1.4683e+10   2.90974e+09  11.97      0.3629                  
2.57841e+10   3.59833e+09  12.16      0.3087                  
4.16434e+10   5.83513e+09  12.43      0.2122                  
6.67911e+10   9.34632e+09  12.61      0.2369                  
1.07375e+11   1.52612e+10  13.26      0.1444                  
1.73309e+11   2.39779e+10  13.34      0.2277                  
2.77191e+11   3.80254e+10  nan      nan                   
4.48464e+11   5.69671e+10  12.9       1.273
];%r<19.4 centrals
% xbrd=[5e9,1e10; 1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11; 3.7e11, 6e11];
data2=[%nan nan nan nan
    1.51429e+10   2.91794e+09  12.24      1.216
    2.61354e+10   3.59945e+09  12.15      1.12                    
4.24234e+10   5.81214e+09  12.76      0.3252                  
6.77309e+10   9.51513e+09  13.1       0.227                   
1.08512e+11   1.54176e+10  13.63      0.1394                  
1.73428e+11   2.33892e+10  13.58      0.2571                       
2.7515e+11   3.71434e+10  14.09      0.2438
nan nan nan nan
]; %N>2
data5=[nan nan nan nan
    1.52189e+10   2.77664e+09  12.92      1.244                   
2.64628e+10   3.86419e+09  12.2       1.607                   
4.41984e+10   5.89932e+09  13.11      0.4577                  
6.90324e+10   9.4982e+09  13.44      0.2624                  
1.09594e+11   1.60071e+10  13.88      0.1412                  
1.74127e+11   2.34103e+10  13.99      0.1942                  
2.73782e+11   3.49934e+10  14.36      0.2142                  
% 4.33703e+11   4.69153e+10  9.908      63.45 
]; %N>5;
dataZ02=[%7.28742e+09   1.45184e+09  11.63      0.5295                  
1.43904e+10   2.89589e+09  11.72      0.5876                  
2.52734e+10   3.55315e+09  12.16      0.3631                  
4.10653e+10   5.77329e+09  12.47      0.2683                  
6.55599e+10   9.16918e+09  12.79      0.2516                  
1.05014e+11   1.43846e+10  13.65      0.1315                  
1.67394e+11   2.12134e+10  13.6       0.2977                  
2.74124e+11   2.95823e+10  13.33      1.34                    
5.09854e+11   7.66566e+10  8.618      nan
]; %z<0.2
dataZ03=[7.36267e+09   1.45174e+09  11.5       0.6058                  
1.46412e+10   2.90471e+09  11.98      0.3591                  
2.57035e+10   3.59548e+09  12.18      0.3014                  
4.13775e+10   5.77911e+09  12.36      0.2411                  
6.58673e+10   9.19272e+09  12.7       0.2228                  
1.05256e+11   1.46938e+10  13.37      0.1417                  
1.69788e+11   2.2613e+10  13.56      0.2031                  
2.7356e+11   3.33571e+10  12.64      2.822                   
4.64753e+11   6.67674e+10  13.35      1.255];%z<0.3
xbrd=[1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11;3.7e11, 6e11];
% xbrd=log10(xbrd);
myfigure;
% sigma_common=0.1;
% for i=1:numel(ms)
%     zwang(i)=medianMh(log10(ms(i)),sigma_common,Pwang4);
%     zling(i)=medianMh(log10(ms(i)),sigma_common,Pling);
%     zmosterZ0(i)=medianMh(log10(ms(i)),sigma_common,PmosterZ0);
%     zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo); 
% end
% zbound=minmax([zwang',zling',zmosterZ0',zguo']);
% h4=area(ms,[zbound(:,1),zbound(:,2)-zbound(:,1)]);hold on;
% set(h4(2),'facecolor','b','edgecolor','w');
% set(h4(1),'facecolor','w','edgecolor','w');
sigma_common=0.2;
for i=1:numel(ms)
    zwang(i)=medianMh(log10(ms(i)),sigma_common,Pwang4);
    zling(i)=medianMh(log10(ms(i)),sigma_common,Pling);
    zmosterZ0(i)=medianMh(log10(ms(i)),sigma_common,PmosterZ0);
    zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo); 
end
zbound=minmax([zwang',zling',zmosterZ0',zguo']);
h3=area(ms,[zbound(:,1),zbound(:,2)-zbound(:,1)]);hold on;
set(h3(2),'facecolor','y','edgecolor','w');
set(h3(1),'facecolor','w','edgecolor','w');
% sigma_common=0.3;
% for i=1:numel(ms)
%     zwang(i)=medianMh(log10(ms(i)),sigma_common,Pwang4);
%     zling(i)=medianMh(log10(ms(i)),sigma_common,Pling);
%     zmosterZ0(i)=medianMh(log10(ms(i)),sigma_common,PmosterZ0);
%     zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo); 
% end
% zbound=minmax([zwang',zling',zmosterZ0',zguo']);
% h5=area(ms,[zbound(:,1),zbound(:,2)-zbound(:,1)]);hold on;
% set(h5(2),'facecolor','c','edgecolor','w');
% set(h5(1),'facecolor','w','edgecolor','w');
set(gca,'layer','top');
% f=grpmock.Mult>5;
% x=logspace(9.8,12,15);
% [xmed,ymed]=skeleton(grpmock.MstarIter(f),log10(grpmock.MIter(f)),x,0.68);
% h55=plot(xmed,ymed,'-','linewidth',2,'color','g');hold on;
f=grpmock.Mult>2;
x=logspace(9,12,15);
[xmed,ymed]=skeleton(grpmock.MstarIter(f),log10(grpmock.MIter(f)),x,0.68);
h22=plot(xmed,ymed,'-','linewidth',2,'color','r');hold on;
f=galmock.CentralSampleIter>0&galmock.Z<0.2;
[xmed,ymed,yl,xm,ym,yerr]=skeleton(galmock.Mstar(f),log10(galmock.Mhalo(f)),x,0.68);
h11=plot(xmed,ymed,'k--');
% h5=ploterr(data5(:,1),data5(:,3),[],data5(:,4),'g^','logx');hold on;
h2=ploterr(data2(:,1),data2(:,3)-0.0,{xbrd(:,1),xbrd(:,2)},data2(:,4),'ro','logx');hold on;
% h1=ploterr(data(:,1),data(:,3)-0.1,{xbrd(:,1),xbrd(:,2)},data(:,4),'kd','logx');hold on;
h1=ploterr(dataZ02(:,1),dataZ02(:,3)-0.,{xbrd(:,1),xbrd(:,2)},dataZ02(:,4),'ks','logx');hold on;set(h1(1),'markersize',15);
% h1=ploterr(dataZ03(:,1),dataZ03(:,3)-0.1,{xbrd(:,1),xbrd(:,2)},dataZ03(:,4),'cs','logx');hold on;
% f=galmock.CentralSampleIter>0&galmock.Z<0.2;
% [xmed,ymed,yl,xm,ym,yerr]=skeleton(galmock.Mstar(f),log10(galmock.Mhalo(f)),x,0.68);
% h11=plot(xmed,ymed,'c--');
% plot(x,12.9-0.1+1.009*(log10(x)-11),'r-');
% plot(x,12.76-0.1+0.84*(log10(x)-11),'r-');
% plot(x,13.06-0.1+1.5*(log10(x)-11),'r-');
xlim([9e9,1e12]);
% xlim([4.5e9,1e12]);
ylim([10.8,15.5]);
xscale('log');
xlabel('$M_\star$[M$_\odot/h^2$]');
ylabel('$\log(M_h$[M$_\odot/h$])');
l=legend([h2(1),h1(1),h22,h11,h3(2)],'$N>2$','$N\geq 1, z<0.2$','Mock $N>2$','Mock $N\geq 1, z<0.2$','Models $\sigma=0.2$');
% h4=plot(ms,ywang4,'c--');
% h5=plot(ms,yling,'m--');
% h6=plot(ms,ymosterZ0,'g--');
% mh=logspace(10,16,10);
% h7=plot(halo2starP(mh,Pguo),log10(mh),'b--'); 

set(l,'location','northwest','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/M-Mstar-central.eps');
%% compare with lensing measurements
dataZ02=[%7.28742e+09   1.45184e+09  11.63      0.5295                  
1.43904e+10   2.89589e+09  11.72      0.5876                  
2.52734e+10   3.55315e+09  12.16      0.3631                  
4.10653e+10   5.77329e+09  12.47      0.2683                  
6.55599e+10   9.16918e+09  12.79      0.2516                  
1.05014e+11   1.43846e+10  13.65      0.1315                  
1.67394e+11   2.12134e+10  13.6       0.2977                  
2.74124e+11   2.95823e+10  13.33      1.34                    
5.09854e+11   7.66566e+10  8.618      nan
]; %z<0.2
xbrd=[1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11;3.7e11, 6e11];
% xbrd=log10(xbrd);
myfigure;% set(gcf,'defaultlinemarkersize',10)
% sigma_common=0.2;
% for i=1:numel(ms)
%     zwang(i)=medianMh(log10(ms(i)),sigma_common,Pwang4);
%     zling(i)=medianMh(log10(ms(i)),sigma_common,Pling);
%     zmosterZ0(i)=medianMh(log10(ms(i)),sigma_common,PmosterZ0);
%     zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo); 
% end
% zbound=minmax([zwang',zling',zmosterZ0',zguo']);
% h3=area(ms,[zbound(:,1),zbound(:,2)-zbound(:,1)]);hold on;
% set(h3(2),'facecolor','y','edgecolor','w');
% set(h3(1),'facecolor','w','edgecolor','w');
% set(gca,'layer','top');

% h0=ploterr(data(:,1),data(:,3)-0.1,{xbrd(:,1),xbrd(:,2)},data(:,4),'rd','logx');hold on; %flux-limited centrals
% h11=ploterr(dataZ02(:,1),dataZ02(:,3)-0.2,[],dataZ02(:,4),'ks','logx');hold on;
% set(h11,'color',[0.8,0.8,0.8]);
% h1=ploterr(dataZ03(:,1),dataZ03(:,3)-0.1,{xbrd(:,1),xbrd(:,2)},dataZ03(:,4),'cs','logx');hold on;
% f=galmock.CentralSampleIter>0&galmock.Z<0.2;
% [xmed,ymed,yl,xm,ym,yerr]=skeleton(galmock.Mstar(f),log10(galmock.Mhalo(f)),x,0.68);
% h11=plot(xmed,ymed,'c--');
% plot(x,12.9-0.1+1.009*(log10(x)-11),'r-');
% plot(x,12.76-0.1+0.84*(log10(x)-11),'r-');
% plot(x,13.06-0.1+1.5*(log10(x)-11),'r-');
% plot(x,log10(1.3e13)+1.49*log10(x/2e11),'r-');%cfht red
% plot(x,log10(0.88e13)+0.83*log10(x/2e11),'b-');%cfht blue
SDSSearly=[0.76,3.18, -3.14, 9.45; %Mstar/1e10Msun, Mh/1e11Msunh, deltaMh
    1.5 4.20  -3.67 6.63
    3.0 4.9  -3.2 4.7
    5.8 14.1 -5.3 5.6 
    11.2 34 -9 10
    21.3 158 -33 37
    39.6 716 -190 123];
SDSSlate=[0.73, 0.02 -0.018 1.56
    1.5 6.6 -4 5.1
    2.9 6.1 -3.4 4.7
    5.6 14 -7 8
    10.8 13 -9 12
    20.5 34 -28 33
    40.4 180 -173 532];
CFHTblue=[0.18, 1.28, -0.33, 0.41 %Mstar/1e10Msun, M200c/1e11Msun, deltaMh
   0.54 2.0 -0.62 0.64
   1.59 9.14 -1.88 2.37
   4.27 26.8 -10.3 11.0];
CFHTred=[0.24 0.03 -0.02 1.90
    0.66 5.68 -1.84 2.16
    1.97 5.81 -1.20 1.67
    5.64 26.3 -2.88 3.23
    13.0 81.2 -8.91 12.1
    22.6 160 -24.2 28.3
    38.6 388 -67.1 90.7
    62.7 174 -167 353];

d=SDSSearly; d(:,1)=d(:,1)*1e10*0.7^2; d(:,2:4)=d(:,2:4)*1e11;
h3=ploterr(d(:,1), log10(d(:,2)), [], {log10(d(:,2)+d(:,3)), log10(d(:,2)+d(:,4))}, 'mo' ,'logx');%set(h3(1), 'markersize',7);
hold on;
d=SDSSlate; d(:,1)=d(:,1)*1e10*0.7^2; d(:,2:4)=d(:,2:4)*1e11;
h2=ploterr(d(:,1), log10(d(:,2)), [], {log10(d(:,2)+d(:,3)), log10(d(:,2)+d(:,4))}, 'gd', 'logx');set(h2(1), 'markersize',8); %, 'markerfacecolor','g');
d=CFHTblue; d(:,1)=d(:,1)*1e10*0.7^2; d(:,2:4)=d(:,2:4)*1e11*0.7/MvcInMvb;
h4=ploterr(d(:,1), log10(d(:,2)), [], {log10(d(:,2)+d(:,3)), log10(d(:,2)+d(:,4))}, 'bv' ,'logx');
d=CFHTred; d(:,1)=d(:,1)*1e10*0.7^2; d(:,2:4)=d(:,2:4)*1e11*0.7/MvcInMvb;
h5=ploterr(d(:,1), log10(d(:,2)), [], {log10(d(:,2)+d(:,3)), log10(d(:,2)+d(:,4))}, 'r^', 'logx');
% h4=plot(ms,ywang4,'c--');
% h5=plot(ms,yling,'m--');
% h6=plot(ms,ymosterZ0,'g--');
% mh=logspace(10,16,10);
% h7=plot(halo2starP(mh,Pguo),log10(mh),'b--'); 
h6=plot(ms,ycosmos,'k--');
h7=plot(ms, yhudson, 'k:');
h1=ploterr(dataZ02(:,1),dataZ02(:,3)-0.0,[],dataZ02(:,4),'ks','logx');hold on;
set(h1,'markersize',18);
xlim([9e9,1e12]);
% xlim([4.5e9,1e12]);
ylim([10.8,15.5]);
xscale('log');
xlabel('$M_\star$[M$_\odot/h^2$]');
ylabel('$\log(M_h$[M$_\odot/h$])');
l=legend([h1(1),h2(1),h3(1),h4(1),h5(1),h7,h6],'GAMA(this work)','SDSS late','SDSS early','CFHT blue(V13)','CFHT red(V13)','CFHT all(H13)','COSMOS');
% h4=plot(ms,ywang4,'c--');
% h5=plot(ms,yling,'m--');
% h6=plot(ms,ymosterZ0,'g--');
% mh=logspace(10,16,10);
% h7=plot(halo2starP(mh,Pguo),log10(mh),'b--'); 

set(l,'location','northwest','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/M-Mstar-central-LensCmp.eps');
%% M=A(M*) comparison, SSFR dep             
data=[7.3685e+09   1.45222e+09  11.51      0.5988                  
1.4683e+10   2.90974e+09  11.97      0.3629                  
2.57841e+10   3.59833e+09  12.16      0.3087                  
4.16434e+10   5.83513e+09  12.43      0.2122                  
6.67911e+10   9.34632e+09  12.61      0.2369                  
1.07375e+11   1.52612e+10  13.26      0.1444                  
1.73309e+11   2.39779e+10  13.34      0.2277                  
2.77191e+11   3.80254e+10  nan      nan                   
4.48464e+11   5.69671e+10  12.9       1.273
];%r<19.4 centrals
xbrd=[5e9,1e10; 1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11; 3.7e11, 6e11];
xbrd=log10(xbrd);
datahigh=[7.29997e+09   1.43735e+09  12.01      0.4588                  
1.442e+10   2.88151e+09  12.51      0.3106                  
2.52676e+10   3.51158e+09  12.8       0.3169                  
3.97499e+10   5.5919e+09  13.11      0.3219                  
6.31942e+10   8.748e+09  5.297     nan                   %266 halos
1.00395e+11   1.35373e+10  13         nan                   %41 halos
1.5887e+11   1.47669e+10  14.22      1.096                   %3 halos in total
2.61078e+11   2.69671e+10  10.28      nan                    %3 halos in total
nan nan nan nan
]; %high SSFR, and restricted to r<19.4
% datahigh=[7.32997e+09   1.44329e+09  11.61      0.6861                  
% 1.44613e+10   2.87489e+09  12.45      0.2842                  
% 2.52509e+10   3.52084e+09  12.7       0.2742                  
% 4.00896e+10   5.66741e+09  12.66      0.3878                  
% 6.34077e+10   8.86156e+09  5.159      nan                   %386 halos
% 1.0049e+11   1.37032e+10  13.39      0.7421                %67 halos  
% 1.76798e+11   2.41545e+10  14.47      0.4019      %15 halos            
% 2.64844e+11   2.68815e+10  10.32      nan %8 halos
% nan nan nan nan]; %high SSFR, FluxHalpha cut
% datahigh=[7.32997e+09   1.44329e+09  11.58      0.7309                  
% 1.44613e+10   2.87489e+09  12.42      0.2989                  
% 2.52509e+10   3.52084e+09  12.7       0.2768                  
% 4.00896e+10   5.66741e+09  12.71      0.3651                  
% 6.34077e+10   8.86156e+09  5.112      nan %667.7                   
% 1.0049e+11   1.37032e+10  13.42      0.7007                  
% 1.76798e+11   2.41545e+10  14.44      0.4014                  
% 2.64844e+11   2.68815e+10  10.43      nan %21.98   
% nan nan nan nan]; %high SSFR, SNCont cut
% datahigh=[7.32997e+09   1.44329e+09  11.62      0.6854                  
% 1.44613e+10   2.87489e+09  12.45      0.2815                  
% 2.52509e+10   3.52084e+09  12.69      0.2757                  
% 4.00896e+10   5.66741e+09  12.68      0.375                   
% 6.34077e+10   8.86156e+09  5.112      nan %667.7                   
% 1.0049e+11   1.37032e+10  13.42      0.7007                  
% 1.76798e+11   2.41545e+10  14.44      0.4014                  
% 2.64844e+11   2.68815e+10  10.43      nan %21.98 
% nan nan nan nan]; %high SSFR, SNComb cut
% datahigh=[7.41597e+09   1.42933e+09  11.87      0.5842                  
% 1.45578e+10   2.8713e+09  12.43      0.3164                  
% 2.54351e+10   3.56045e+09  12.69      0.2761                  
% 4.07431e+10   5.75681e+09  12.14      0.665                   
% 6.50974e+10   9.13563e+09  5.072      nan %459                     
% 1.0454e+11   1.49481e+10  13.5       0.3178                  
% 1.71442e+11   2.26322e+10  14.06      0.2977                  
% 2.67164e+11   2.68209e+10  11.36     nan % 6.008                   
% 4.90697e+11   6.3415e+10  14.11      0.687      
% ]; %high LumHalpha
datacomp=[7.44261e+09   1.46453e+09  8.88       nan                   
1.48636e+10   2.91526e+09  11.35      1.267                   
2.59596e+10   3.61056e+09  11.96      0.4464                  
4.18871e+10   5.82114e+09  12.35      0.2438                  
6.6924e+10   9.3412e+09  12.64      0.23                    
1.07442e+11   1.52613e+10  13.25      0.1453                  
1.73332e+11   2.3983e+10  13.34      0.2288                  
2.77271e+11   3.80553e+10  5.082     nan                   
4.48464e+11   5.69671e+10  12.9       1.273 ]; %complementary to high SSFR, within r<19.4 centrals
datazbound=[7.36217e+09   1.44739e+09  11.53      0.5673                  
1.46087e+10   2.90045e+09  11.85      0.4284                  
2.5687e+10   3.5878e+09  12.21      0.2809                  
4.13266e+10   5.79956e+09  12.38      0.2283                  
6.58923e+10   9.15656e+09  12.63      0.2415                  
1.05298e+11   1.46054e+10  13.36      0.1415                  
1.69995e+11   2.26231e+10  13.54      0.2051                  
2.7353e+11   3.49233e+10  12.65      2.08                    
4.63165e+11   6.3437e+10  13.37      0.9454]; %the same redshift range as the ssfr subsample
databluehigh=[7.2836e+09   1.43616e+09  11.91      0.5504                  
1.42952e+10   2.86541e+09  12.59      0.3139                  
2.49872e+10   3.41752e+09  12.8       0.3918                  
3.94474e+10   5.47259e+09  13.27      0.3479                  
6.31573e+10   8.71878e+09  13.44      0.7397                  
9.52068e+10   9.99711e+09  8.435      nan                   
1.64675e+11   1.50332e+10  14.27      1.03                    
2.78027e+11   1.51316e+10  10.34      nan
nan nan nan nan]; %blue and high
datablueOrhigh=[7.27395e+09   1.43953e+09  11.48      0.7363                  
1.44586e+10   2.89771e+09  12.01      0.4935                  
2.55198e+10   3.57771e+09  12.13      0.6385                  
4.07444e+10   5.77469e+09  12.55      0.4249                  
6.58821e+10   9.35674e+09  5.002      nan                   
1.08407e+11   1.5415e+10  10.2       nan                   
1.7403e+11   2.49996e+10  13.63      0.4905                  
2.80344e+11   3.72752e+10  6.914      nan                   
4.6118e+11   6.19565e+10  11.25      nan]; %blue or high
datahighz=[7.30049e+09   1.45166e+09  11.96      0.6161                  
1.46217e+10   2.87625e+09  12.49      0.3783                  
2.59301e+10   3.58184e+09  12.23      0.7441                  
4.15344e+10   5.82956e+09  12.64      0.3616                  
6.67339e+10   9.31684e+09  12.69      0.572                   
1.07675e+11   1.54596e+10  12.46      1.126                   
1.74241e+11   2.41788e+10  13.53      0.4475                  
2.76378e+11   3.8721e+10  8.453      nan                  
4.39403e+11   5.11245e+10  13.25      1.325 ];
dataactive=[7.32956e+09   1.43946e+09  11.91      0.5698                  
1.44526e+10   2.87549e+09  12.74      0.2417                  
2.52757e+10   3.51309e+09  12.75      0.3466                  
3.97927e+10   5.59366e+09  13.13      0.3157                  
6.3215e+10   8.75789e+09  5.294      nan %829.6                   
9.98413e+10   1.32395e+10  8.854      nan %88.25                   
1.63484e+11   1.62241e+10  14.98      0.595                   
2.61078e+11   2.69671e+10  10.28      nan %25.58
nan nan nan nan];
datapassive=[7.36897e+09   1.44751e+09  11.33      0.7568                  
1.51124e+10   2.8567e+09  11.77      0.5317                  
2.6223e+10   3.59262e+09  12.17      0.3295                  
4.16283e+10   5.68136e+09  12.53      0.1983                  
6.59916e+10   9.22755e+09  12.66      0.228                   
1.0261e+11   1.32545e+10  13.32      0.1395                  
1.66816e+11   2.24307e+10  13.22      0.2761                  
2.76074e+11   2.46063e+10  5.191      873.1                   
4.8631e+11   6.01643e+10  12.37      2.122 ];
dataactive1R=[7.32956e+09   1.43946e+09  12.08      0.536                   
1.44526e+10   2.87549e+09  12.66      0.3163                  
2.52757e+10   3.51309e+09  12.76      0.3976                  
3.97927e+10   5.59366e+09  13.14      0.3733                  
6.3215e+10   8.75789e+09  5.162      nan %934.5                   
9.98413e+10   1.32395e+10  13.94      1.069                   
1.63484e+11   1.62241e+10  8.742      nan %111.4                   
2.61078e+11   2.69671e+10  10.47      nan %21.75
nan nan nan nan];
dataactive1R03offset=[7.32956e+09   1.43946e+09  12.09      0.5292                  
1.44526e+10   2.87549e+09  12.66      0.3163                  
2.52757e+10   3.51309e+09  12.74      0.4097                  
3.97927e+10   5.59366e+09  13.17      0.3579                  
6.3215e+10   8.75789e+09  5.162      nan %934.5                   
9.98413e+10   1.32395e+10  13.94      1.069                   
1.63484e+11   1.62241e+10  8.742      nan %111.4                   
2.61078e+11   2.69671e+10  9.082      nan %55.84 
nan nan nan nan]; %offsetcut=0.3
dataactiveNewEstimate=[7.32956e+09   1.43946e+09  12.02      0.567                   
1.44526e+10   2.87549e+09  12.74      0.2866                  
2.52757e+10   3.51309e+09  12.65      0.4323                  
3.97927e+10   5.59366e+09  13.11      0.3649                  
6.3215e+10   8.75789e+09  5.121     nan % 996.6                   
9.98413e+10   1.32395e+10  13.45     nan % 1.856                   
1.63484e+11   1.62241e+10  6.498      nan %1074                    
2.61078e+11   2.69671e+10  10.26      nan %26.15
nan nan nan nan]; %using the active M-M* relation to estimate LumMass and Rvir, 2Rv/0offset cut.
myfigure;
sigma_common=0.2;
for i=1:numel(ms)
    zwang(i)=medianMh(log10(ms(i)),sigma_common,Pwang4);
    zling(i)=medianMh(log10(ms(i)),sigma_common,Pling);
    zmosterZ0(i)=medianMh(log10(ms(i)),sigma_common,PmosterZ0);
    zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo); 
end
zbound=minmax([zwang',zling',zmosterZ0',zguo']);
h3=area(log10(ms),[zbound(:,1),zbound(:,2)-zbound(:,1)]);hold on;
set(h3(2),'facecolor','y','edgecolor','w');
set(h3(1),'facecolor','w','edgecolor','w');
set(gca,'layer','top');
%----------mock predictions----------------
x=linspace(8,12,12);
%---raw mocks---- like shit
% f=rawmock.z<0.5;
% [xmed,ymed]=skeleton(rawmock.sm(f),log10(galmock.Mhalo(f)),x,0.68);
% h2=plot(xmed,ymed,'-d','linewidth',2,'color','r');hold on;
% f=rawmock.z<0.5&rawmock.SFR./rawmock.sm>10^-1.5;
% [xmed,ymed]=skeleton(rawmock.sm(f),log10(galmock.Mhalo(f)),x,0.68);
% h2=plot(xmed,ymed,'-d','linewidth',2,'color','b');hold on;
%---real centrals-----
% f=galmock.is_central>0;
% x=logspace(8,12,15);
% [xmed,ymed,yl,xm,ym,yerr]=skeleton(galmock.Mstar(f),log10(galmock.Mhalo(f)),x,0.68);
% h22=plot(xmed,ymed,'--','linewidth',2,'color','r');hold on;
% for vol=1:9
%     load(['G3Cv4up8/mockcat_',num2str(vol),'.mat']);
%-------standard mock--------------
f=galmock.CentralSampleIter>0; %about 90% are real centrals
[xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(galmock.Mstar(f)),log10(galmock.Mhalo(f)),x,0.68);
h2=plot(xmed,ymed,'-','linewidth',2,'color','r');hold on;
% end
% plot(xmed,yl,'r--');
f=galmock.CentralSampleIter>0&galmock.SFR./galmock.Mstar>10^-1.5;
[xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(galmock.Mstar(f)),log10(galmock.Mhalo(f)),x,0.68);
h22=plot(xmed,ymed,'-','linewidth',2,'color','b');hold on;
% plot(xmed,yl,'b--');
% f=(galmock.EnvLevel==0|galmock.IsIterCen>0)&galmock.SFR./galmock.Mstar<10^-1.5;
% x=logspace(8,12,15);
% [xmed,ymed,yl,xm,ym,yerr]=skeleton(galmock.Mstar(f),log10(galmock.Mhalo(f)),x,0.68);
% h222=plot(xmed,ymed,'-','linewidth',2,'color','g');hold on;
%---------Eagle----------------
% f=Eagle.SSFR<10.^-2; 
% [xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(Eagle.Mgal(f)),log10(Eagle.Mhalo(f)),x,0.68);
% h21=plot(xmed,ymed,'-','linewidth',2,'color','m');hold on;
% plot(xmed,yl,'r--');
% f=Eagle.SSFR>10^-2;
% [xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(Eagle.Mgal(f)),log10(Eagle.Mhalo(f)),x,0.68);
% h221=plot(xmed,ymed,'-','linewidth',2,'color','g');hold on;
% plot(xmed,yl,'b--');
% ---------adding noise in Mstar------------------
galmock.Mstar2=galmock.Mstar.*10.^(normrnd(0.05,0.4,size(galmock.Mstar)));
galmock.SFR2=galmock.SFR.*10.^(normrnd(0.12,0.37,size(galmock.Mstar)));
% % galmock.SFR2=galmock.SFR2(randperm(numel(galmock.Mstar))); %permute to remove SFR dependence
f=galmock.CentralSampleIter>0; %about 90% are real centrals
[xmed,ymed]=skeleton(log10(galmock.Mstar2(f)),log10(galmock.Mhalo(f)),x,0.68);
h20=plot(xmed,ymed,'--','linewidth',2,'color','r');hold on;
f=galmock.CentralSampleIter>0&galmock.SFR2./galmock.Mstar2>10^-1.5;
[xmed,ymed]=skeleton(log10(galmock.Mstar2(f)),log10(galmock.Mhalo(f)),x,0.68);
h220=plot(xmed,ymed,'--','linewidth',2,'color','b');hold on;
% ------------volume limited mock----------------
% magcut=19.4; %almost no effect in applying a flux limit or not, in the mocks.
% tmpmock=bower;
% f=tmpmock.rmag<magcut;
% [xmed,ymed]=skeleton(tmpmock.sm(f),log10(tmpmock.mh(f,1)),x,0.68);
% h220=plot(xmed,ymed,'x--','linewidth',2,'color','r');hold on;
% f=tmpmock.rmag<magcut&tmpmock.sfr./tmpmock.sm>10^-1.5&tmpmock.z<0.31;
% [xmed,ymed]=skeleton(tmpmock.sm(f),log10(tmpmock.mh(f,1)),x,0.68);
% h220=plot(xmed,ymed,'x--','linewidth',2,'color','b');hold on;
% tmpmock=guo;
% f=tmpmock.rmag<magcut;
% [xmed,ymed]=skeleton(tmpmock.sm(f),log10(tmpmock.mh(f,1)),x,0.68);
% h220=plot(xmed,ymed,'--','linewidth',2,'color','r');hold on;
% f=tmpmock.rmag<magcut&tmpmock.sfr./tmpmock.sm>10^-1.5&tmpmock.z<0.31;
% [xmed,ymed]=skeleton(tmpmock.sm(f),log10(tmpmock.mh(f,1)),x,0.68);
% h220=plot(xmed,ymed,'--','linewidth',2,'color','b');hold on;
% tmpmock=font;
% f=tmpmock.rmag<magcut;
% [xmed,ymed]=skeleton(tmpmock.sm(f),log10(tmpmock.mh(f,1)),x,0.68);
% h220=plot(xmed,ymed,'<--','linewidth',2,'color','r');hold on;
% f=tmpmock.rmag<magcut&tmpmock.sfr./tmpmock.sm>10^-1.5&tmpmock.z<0.31;
% [xmed,ymed]=skeleton(tmpmock.sm(f),log10(tmpmock.mh(f,1)),x,0.68);
% h220=plot(xmed,ymed,'<--','linewidth',2,'color','b');hold on;
%---------HOD models---------------------
% plot(xmed,13.55+0.81*log10(xmed/1e12),'k-');
% plot(xmed,13.85+0.92*log10(xmed/1e12),'k--');
% h4=plot(ms,ywang4,'c--');
% h5=plot(ms,yling,'m--');
% h6=plot(ms,ymosterZ0,'g--');
% mh=logspace(10,16,10);
% h7=plot(log10(halo2starP(mh,Pguo)),log10(mh),'b--'); 


h1=ploterr(log10(data(:,1)),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro');hold on;
% h11=ploterr(datacomp(:,1),datacomp(:,3),{xbrd(:,1),xbrd(:,2)},datacomp(:,4),'gd','logx');hold on;
% h11=ploterr(log10(datacomp(:,1)),datacomp(:,3),[],datacomp(:,4),'gd');hold on;
% h11=ploterr(datazbound(:,1),datazbound(:,3),{xbrd(:,1),xbrd(:,2)},datazbound(:,4),'gx','logx');hold on;
f=datahigh(:,1)<1e12;
h111=ploterr(log10(datahigh(f,1)),datahigh(f,3),{xbrd(f,1),xbrd(f,2)},datahigh(f,4),'bs','logx');hold on;
% h111=ploterr(log10(dataactive(f,1)),dataactive(f,3),{xbrd(f,1),xbrd(f,2)},dataactive(f,4),'mo','logx');hold on;
% h=ploterr(log10(datahighz(:,1)),datahighz(:,3),[],[],'bx');set(h(1),'markersize',15);
% plot(x,log10(0.88e13)+0.83*(x-log10(2e11)),'b-'); %cfht fit, blue
% plot(x,log10(1.3e13)+1.49*(x-log10(2e11)),'r-'); %cfht fit, red
% plot(x,15.2+1.58*(x-12),'b-');
% plot(x,15.01+1.43*(x-12),'r-');
% plot(x,13.96+0.93*(x-12),'g-');
% plot(x,13.75+0.78*(x-12),'g--');
% plot(x,13.59+0.85*(x-12),'r-');
xlim([9.5,12]);
ylim([10,15.5]);
% xscale('log');
xlabel('$\log(M_\star$[M$_\odot/h^2$])');
ylabel('$\log(M_h$[M$_\odot/h$])');
l=legend([h1(1),h111(1),h11(1),h2,h22,h3(2)],'All Central','Active Central','Inactive Central','Mock Central','Mock Active','$\sigma=0.2$ HOD');
% l=legend([h1(1),h2,h22,h222,h4,h5,h6,h7,h3(2)],'Data','Mock','Mock Active','Mock Passive','WangL13','WangLY13','Moster13','Guo10','All $\sigma=0.2$');
set(l,'location','northwest','interpreter','latex');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/M-Mstar-central-randomnoise.eps');
%% datacut dependence         
data=[7.3685e+09   1.45222e+09  11.51      0.5988                  
1.4683e+10   2.90974e+09  11.97      0.3629                  
2.57841e+10   3.59833e+09  12.16      0.3087                  
4.16434e+10   5.83513e+09  12.43      0.2122                  
6.67911e+10   9.34632e+09  12.61      0.2369                  
1.07375e+11   1.52612e+10  13.26      0.1444                  
1.73309e+11   2.39779e+10  13.34      0.2277                  
2.77191e+11   3.80254e+10  nan      nan                   
4.48464e+11   5.69671e+10  12.9       1.273
];%r<19.4 centrals
xbrd=[5e9,1e10; 1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11; 3.7e11, 6e11];
xbrd=log10(xbrd);
datahigh=[7.29997e+09   1.43735e+09  12.01      0.4588                  
1.442e+10   2.88151e+09  12.51      0.3106                  
2.52676e+10   3.51158e+09  12.8       0.3169                  
3.97499e+10   5.5919e+09  13.11      0.3219                  
6.31942e+10   8.748e+09  5.297     nan                   %266 halos
1.00395e+11   1.35373e+10  13         nan                   %41 halos
1.5887e+11   1.47669e+10  14.22      1.096                   %3 halos in total
2.61078e+11   2.69671e+10  10.28      nan                    %3 halos in total
nan nan nan nan
]; %high SSFR, and restricted to r<19.4
datacomp=[7.44261e+09   1.46453e+09  8.88       nan                   
1.48636e+10   2.91526e+09  11.35      1.267                   
2.59596e+10   3.61056e+09  11.96      0.4464                  
4.18871e+10   5.82114e+09  12.35      0.2438                  
6.6924e+10   9.3412e+09  12.64      0.23                    
1.07442e+11   1.52613e+10  13.25      0.1453                  
1.73332e+11   2.3983e+10  13.34      0.2288                  
2.77271e+11   3.80553e+10  5.082     nan                   
4.48464e+11   5.69671e+10  12.9       1.273 ]; %complementary to high SSFR, within r<19.4 centrals
datahighall=[7.44261e+09   1.46453e+09  11.62      0.5988                  
1.48636e+10   2.91526e+09  12.35      0.2955                  
2.59596e+10   3.61056e+09  12.42      0.3519                  
4.18871e+10   5.82114e+09  12.6       0.3539                  
6.6924e+10   9.3412e+09  12.33      1.1                     
1.07442e+11   1.52613e+10  13.61      0.3005                  
1.73332e+11   2.3983e+10  8.545      nan %60.87                   
2.77271e+11   3.80553e+10  10.4       nan %24.52                   
4.48464e+11   5.69671e+10  8.93       nan ];%94.37 ] % no FlagSFR cut
datalowall=[7.44261e+09   1.46453e+09  11.56      1.413                   
1.48636e+10   2.91526e+09  12.04      0.5382                  
2.59596e+10   3.61056e+09  12.17      0.4724                  
4.18871e+10   5.82114e+09  12.07      0.4889                  
6.6924e+10   9.3412e+09  12.72      0.2588                  
1.07442e+11   1.52613e+10  13.23      0.1907                  
1.73332e+11   2.3983e+10  13.19      0.3441                  
2.77271e+11   3.80553e+10  5.355     nan % 954.2                   
4.48464e+11   5.69671e+10  12.97      1.862 ]; % only require SSFR>0, and SSFR<10^-1.5, and central
dataactive=[7.32956e+09   1.43946e+09  11.91      0.5698                  
1.44526e+10   2.87549e+09  12.74      0.2417                  
2.52757e+10   3.51309e+09  12.75      0.3466                  
3.97927e+10   5.59366e+09  13.13      0.3157                  
6.3215e+10   8.75789e+09  5.294      nan %829.6                   
9.98413e+10   1.32395e+10  8.854      nan %88.25                   
1.63484e+11   1.62241e+10  14.98      0.595                   
2.61078e+11   2.69671e+10  10.28      nan %25.58
nan nan nan nan];
% datapassive=[7.36897e+09   1.44751e+09  11.33      0.7568                  
% 1.51124e+10   2.8567e+09  11.77      0.5317                  
% 2.6223e+10   3.59262e+09  12.17      0.3295                  
% 4.16283e+10   5.68136e+09  12.53      0.1983                  
% 6.59916e+10   9.22755e+09  12.66      0.228                   
% 1.0261e+11   1.32545e+10  13.32      0.1395                  
% 1.66816e+11   2.24307e+10  13.22      0.2761                  
% 2.76074e+11   2.46063e+10  5.191      nan %873.1                   
% 4.8631e+11   6.01643e+10  12.37      2.122 ];%wrong, used ActiveSample!=0
datapassive=[7.36897e+09   1.44751e+09  12.77      0.4553                  
1.51124e+10   2.8567e+09  12.65      0.3879                  
2.6223e+10   3.59262e+09  12         0.9827                  
4.16283e+10   5.68136e+09  5.665      nan %346.5                   
6.59916e+10   9.22755e+09  5.148      nan %585.9                   
1.0261e+11   1.32545e+10  12.12      2.166                   
1.66816e+11   2.24307e+10  14.25      0.2838                  
2.76074e+11   2.46063e+10  9.358      nan %47.12                   
4.8631e+11   6.01643e+10  14.33      0.579 
];
datapassive1R=[7.36897e+09   1.44751e+09  11.89      0.442                   
1.51124e+10   2.8567e+09  11.23      1.218                   
2.6223e+10   3.59262e+09  12.17      0.3588                  
4.16283e+10   5.68136e+09  12.66      0.1894                  
6.59916e+10   9.22755e+09  12.85      0.2053                  
1.0261e+11   1.32545e+10  13.35      0.1555                  
1.66816e+11   2.24307e+10  13.12      0.3403                  
2.76074e+11   2.46063e+10  5.209      nan %840                     
4.8631e+11   6.01643e+10  12.73      1.722    ];
datapassive1R03offset=[7.36897e+09   1.44751e+09  11.88      0.4518                  
1.51124e+10   2.8567e+09  11.15      1.414                   
2.6223e+10   3.59262e+09  12.13      0.39                    
4.16283e+10   5.68136e+09  12.69      0.1737                  
6.59916e+10   9.22755e+09  12.87      0.2                     
1.0261e+11   1.32545e+10  13.34      0.1578                  
1.66816e+11   2.24307e+10  13.09      0.3563                  
2.76074e+11   2.46063e+10  5.252      nan %882.3                   
4.8631e+11   6.01643e+10  12.58      2.157];
dataactive1R=[7.32956e+09   1.43946e+09  12.08      0.536                   
1.44526e+10   2.87549e+09  12.66      0.3163                  
2.52757e+10   3.51309e+09  12.76      0.3976                  
3.97927e+10   5.59366e+09  13.14      0.3733                  
6.3215e+10   8.75789e+09  5.162      nan %934.5                   
9.98413e+10   1.32395e+10  13.94      1.069                   
1.63484e+11   1.62241e+10  8.742      nan %111.4                   
2.61078e+11   2.69671e+10  10.47      nan %21.75
nan nan nan nan];
dataactive1R03offset=[7.32956e+09   1.43946e+09  12.09      0.5292                  
1.44526e+10   2.87549e+09  12.66      0.3163                  
2.52757e+10   3.51309e+09  12.74      0.4097                  
3.97927e+10   5.59366e+09  13.17      0.3579                  
6.3215e+10   8.75789e+09  5.162      nan %934.5                   
9.98413e+10   1.32395e+10  13.94      1.069                   
1.63484e+11   1.62241e+10  8.742      nan %111.4                   
2.61078e+11   2.69671e+10  9.082      nan %55.84 
nan nan nan nan]; %offsetcut=0.3
dataactiveNewEstimate=[7.32956e+09   1.43946e+09  12.02      0.567                   
1.44526e+10   2.87549e+09  12.74      0.2866                  
2.52757e+10   3.51309e+09  12.65      0.4323                  
3.97927e+10   5.59366e+09  13.11      0.3649                  
6.3215e+10   8.75789e+09  5.121     nan % 996.6                   
9.98413e+10   1.32395e+10  13.45     nan % 1.856                   
1.63484e+11   1.62241e+10  6.498      nan %1074                    
2.61078e+11   2.69671e+10  10.26      nan %26.15
nan nan nan nan]; %using the active M-M* relation to estimate LumMass and Rvir, 2Rv/0offset cut.
dataFlagSFR=[7.32091e+09   1.44389e+09  12.13      0.3729                  
1.45722e+10   2.89661e+09  12.48      0.257                   
2.56638e+10   3.5593e+09  12.67      0.2523                  
4.10749e+10   5.76057e+09  12.48      0.3939                  
6.55613e+10   9.20045e+09  5.702      nan %386.4                   
1.0186e+11   1.31902e+10  12.73      0.9396                  
1.65392e+11   2.22663e+10  14.22      0.2653                  
2.74587e+11   2.58565e+10  9.641      nan %46.94                   
4.8631e+11   6.01643e+10  14.33      0.579]; %FlagSFR>=1
dataFlagSFRcomp=[7.16252e+09   1.43337e+09  7.082      nan %147.8                   
1.4082e+10   2.87853e+09  6.303      nan %157.7                   
2.51062e+10   3.47152e+09  11.79      0.7015                  
4.01319e+10   5.64045e+09  12.51      0.2667                  
6.35634e+10   8.43211e+09  12.89      0.2279                  
1.01672e+11   1.00588e+10  13.64      0.1324                  
1.66153e+11   2.69242e+10  13.49      0.3539                  
2.77542e+11   1.94806e+10  14.19      0.3298  
nan nan nan nan]; %FlagSFR<1
% dataFlagSFR=[7.32091e+09   1.44389e+09  12.1       0.4267                  
% 1.45722e+10   2.89661e+09  12.72      0.2024                  
% 2.56638e+10   3.5593e+09  12.49      0.3607                  
% 4.10749e+10   5.76057e+09  12.31      0.5904                  
% 6.55613e+10   9.20045e+09  5.195      nan %510.9                   
% 1.0186e+11   1.31902e+10  12.01      2.708                   
% 1.65392e+11   2.22663e+10  14.29      0.2677                  
% 2.74587e+11   2.58565e+10  9.332      nan %67.93                   
% 4.8631e+11   6.01643e+10  14.33      0.579 ]; %FlagSFR>=2
% dataFlagSFRcomp=[7.16252e+09   1.43337e+09  10.31      4.015                   
% 1.4082e+10   2.87853e+09  6.369      nan %144.8                   
% 2.51062e+10   3.47152e+09  12.13      0.4193                  
% 4.01319e+10   5.64045e+09  12.5       0.2678                  
% 6.35634e+10   8.43211e+09  12.87      0.2334                  
% 1.01672e+11   1.00588e+10  13.65      0.1319                  
% 1.66153e+11   2.69242e+10  13.49      0.337                   
% 2.77542e+11   1.94806e+10  14.19      0.3298 
% nan nan nan nan]; %FlagSFR<2
dataZ02high=[7.10085e+09   1.41445e+09  12.28      0.3681                  
1.37583e+10   2.82064e+09  12.65      0.3492                  
2.49631e+10   3.5233e+09  13.02      0.423                   
3.8196e+10   4.80738e+09  4.664      nan %835.3   92groups                   
5.9582e+10   5.99952e+09  9.294      nan %77.65                   
1.04278e+11   1.27987e+10  7.403      nan %202.6                   
1.49642e+11   0  8.105      nan %166.3                   
2.93159e+11   0  10.69      nan %17.59
nan nan nan nan]; %high SSFR, z<0.2, FlagSFR>0
dataZ02active=[7.16252e+09   1.43337e+09  12.25      0.4294                  
1.4082e+10   2.87853e+09  13.04      0.237                   
2.51062e+10   3.47152e+09  12.9       0.5226           %308groups       
4.01319e+10   5.64045e+09  4.678      nan %886.5 , 79groups                  
6.35634e+10   8.43211e+09  9.26       nan %80.37 , 11groups                  
1.01672e+11   1.00588e+10  7.458      nan %165.9                   
1.66153e+11   2.69242e+10  13         1                       
2.77542e+11   1.94806e+10  10.69      nan %17.59 
nan nan nan nan
]; %limited to z<0.2, volume limited for Mstar>1e10, activesample
dataZ02passive=[7.16252e+09   1.43337e+09  12.75      0.4693                  
1.4082e+10   2.87853e+09  12.24      0.7671                  
2.51062e+10   3.47152e+09  11.41      2.602                   
4.01319e+10   5.64045e+09  12.25      2.265                   
6.35634e+10   8.43211e+09  5.947      nan %762.6   , 69groups              
1.01672e+11   1.00588e+10  13.67      0.8964                  
1.66153e+11   2.69242e+10  14.59      0.3776                  
2.77542e+11   1.94806e+10  9.496      nan %52.97
nan nan nan nan
];
dataZ02=[7.28742e+09   1.45184e+09  11.63      0.5295                  
1.43904e+10   2.89589e+09  11.72      0.5876                  
2.52734e+10   3.55315e+09  12.16      0.3631                  
4.10653e+10   5.77329e+09  12.47      0.2683                  
6.55599e+10   9.16918e+09  12.79      0.2516                  
1.05014e+11   1.43846e+10  13.65      0.1315                  
1.67394e+11   2.12134e+10  13.6       0.2977                  
2.74124e+11   2.95823e+10  13.33      1.34                    
5.09854e+11   7.66566e+10  8.618      nan]; %z<0.2
dataZ02activeFluxCut=[7.09709e+09   1.41534e+09  11.88      0.5299                  
1.37095e+10   2.79091e+09  12.52      0.3162                  
2.49795e+10   3.54072e+09  12.73      0.3367                  
3.8551e+10   4.84243e+09  12.01      1.166  %386 groups                   
5.97559e+10   6.23727e+09  13.05      0.7521 %91 groups                 
9.81943e+10   8.38876e+09  10.7      nan %err=21.48
nan nan nan nan];
dataZ02activeCombCut=[7.15494e+09   1.43156e+09  11.89      0.5284                  
1.40422e+10   2.85425e+09  12.52      0.3122                  
2.50575e+10   3.55092e+09  12.72      0.3396                  
4.05202e+10   5.76397e+09  12.1       1.012                   
6.21844e+10   7.95268e+09  11.86      nan%7.407 98 groups                  
9.76211e+10   8.53428e+09  13         nan%4.422                   
1.41429e+11   0  14.35      1.572                   
2.93159e+11   0  10.69      nan%17.59      
nan nan nan nan];
dataZ02activeContCut=[7.15187e+09   1.43135e+09  11.82      0.5829                  
1.40376e+10   2.85092e+09  12.51      0.3195                  
2.50479e+10   3.54699e+09  12.75      0.3302                  
4.05233e+10   5.77098e+09  12.14      0.9644                  
6.21844e+10   7.95268e+09  11.86      nan%7.407 98groups                   
9.76211e+10   8.53428e+09  13         nan%4.422                   
1.41429e+11   0  14.35      1.572                   
2.93159e+11   0  10.69      nan %17.59 
nan nan nan nan];
dataFieldHighZ02=[7.09747e+09   1.41784e+09  12.27      0.3983                  
1.36855e+10   2.83107e+09  12.56      0.4399                  
2.48861e+10   3.47542e+09  12.9       0.6053                  
3.8682e+10   5.01944e+09  5.909      nan %862.4                   
6.0281e+10   6.55558e+09  8.559      nan%114.6                   
9.81943e+10   8.38876e+09  7.458      nan%165.9                   
1.49642e+11   0  8.105      nan%166.3                   
2.93159e+11   0  10.69      nan%17.59    
nan nan nan nan]; %z<0.2 Field centrals (EnvLevel==0 and CenSampleIter>0), high ssfr
dataFieldHigh=[7.30657e+09   1.43933e+09  12.06      0.4566                  
1.44045e+10   2.88908e+09  12.4       0.383                   
2.52052e+10   3.48343e+09  12.73      0.3876                  
3.95086e+10   5.56171e+09  13.24      0.2964                  
6.28997e+10   8.49361e+09  12.69      1.653                   
1.00759e+11   1.35262e+10  8.244      nan %63.93                   
1.64675e+11   1.50332e+10  14.27      1.03                    
2.61078e+11   2.69671e+10  10.28      nan%25.58 
nan nan nan nan]; %Field centrals (EnvLevel==0 and CenSampleIter>0), high ssfr
dataGroupHighZ02=[7.13093e+09   1.38362e+09  12.32      0.9711                  
1.41634e+10   2.72666e+09  12.91      0.5687                  
2.52212e+10   3.66758e+09  13.2       0.5971                  
3.68191e+10   3.82643e+09  4.593      nan %911.4                   
5.7485e+10   3.03102e+09  10.16      nan %28.48                   
1.22528e+11   0  14.99      0.6348       
nan nan nan nan
nan nan nan nan]; %z<0.2 Group Centrals (EnvLevel>0 and CenSampleIter>0), high SSFR
dataGroupHigh=[7.21794e+09   1.4091e+09  11.51      2.296                   
1.45545e+10   2.80457e+09  12.92      0.484                   
2.56843e+10   3.64784e+09  12.99      0.5379                  
4.08173e+10   5.63286e+09  4.601      nan%907.4  216groups                  
6.40029e+10   9.36386e+09  5.822      nan %955.7                   
9.9929e+10   1.3537e+10  14.16      0.6667                  
1.4726e+11   0  8.424      nan%266.2 
nan nan nan nan]; %Group Centrals (EnvLevel>0 and CenSampleIter>0), high SSFR
dataGroup=[7.36856e+09   1.439e+09  10.81      5.043                   
1.49508e+10   2.94332e+09  12.38      0.445                   
2.59873e+10   3.69491e+09  11.78      1.282                   
4.21248e+10   5.78444e+09  12.38      0.3562                  
6.72419e+10   9.53521e+09  12.73      0.2806                  
1.07723e+11   1.50815e+10  13.38      0.1597                  
1.73102e+11   2.32395e+10  13.51      0.2278                  
2.74603e+11   3.73701e+10  13.65      0.4133                  
4.48152e+11   5.58734e+10  13.17      1.574];%All groups
dataField=[7.36951e+09   1.4538e+09  11.55      0.6003                  
1.46399e+10   2.90146e+09  11.76      0.5756                  
2.57438e+10   3.57693e+09  12.25      0.3086                  
4.15308e+10   5.84611e+09  12.45      0.2634                  
6.66129e+10   9.26481e+09  12.43      0.4503                  
1.07199e+11   1.5349e+10  13.01      0.3173                  
1.7329e+11   2.43305e+10  12.83      0.8221                  
2.78659e+11   3.83228e+10  5.59       nan %698.9                   
4.48586e+11   5.73868e+10  12.72      1.905 ]; %Field centrals
dataGroupZ02=[7.33998e+09   1.43439e+09  11.74      1.261                   
1.4886e+10   2.95764e+09  12.29      0.5198                  
2.57734e+10   3.64964e+09  11.83      1.181                   
4.16951e+10   5.86541e+09  12.37      0.4069                  
6.63372e+10   9.3968e+09  12.9       0.2671                  
1.05858e+11   1.46249e+10  13.68      0.14                    
1.67697e+11   2.09796e+10  13.64      0.2942                  
2.72926e+11   3.23325e+10  14.35      0.2843                  
4.33198e+11   0  8.292      113.3    ]; %z<0.2 All groups (N>=2)
dataFieldZ02=[7.2813e+09   1.45424e+09  11.6       0.5851                  
1.42795e+10   2.86977e+09  11.16      1.928                   
2.50816e+10   3.50044e+09  12.26      0.3778                  
4.06888e+10   5.69124e+09  12.57      0.3438                  
6.43378e+10   8.64639e+09  12.53      0.6308                  
1.02034e+11   1.32506e+10  13.48      0.3778                  
1.64312e+11   2.1737e+10  6.131      nan %240                     
2.79942e+11   1.72115e+10  10.01      nan %30.61                   
5.86511e+11   0  10.37      nan %60.35     
]; %z<0.2 Field Centrals
myfigure;

%-------standard mock--------------
f=galmock.CentralSampleIter>0;%&galmock.Z<0.2; %about 90% are real centrals
x=linspace(8,12,12);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(galmock.Mstar(f)),log10(galmock.Mhalo(f)),x,0.68);
h2=plot(xmed,ymed,'-','linewidth',2,'color','r');hold on;
% end
% plot(xmed,yl,'r.-');
f=galmock.CentralSampleIter>0&galmock.SFR./galmock.Mstar>10^-1.5;
[xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(galmock.Mstar(f)),log10(galmock.Mhalo(f)),x,0.68);
h22=plot(xmed,ymed,'-','linewidth',2,'color','b');hold on;

% h1=ploterr(log10(dataZ02(:,1)),dataZ02(:,3),{xbrd(:,1),xbrd(:,2)},dataZ02(:,4),'ko','logx');hold on;
% h1=ploterr(log10(dataFlagSFR(:,1)),dataFlagSFR(:,3),{xbrd(:,1),xbrd(:,2)},dataFlagSFR(:,4),'go','logx');hold on;
% h1=ploterr(log10(data(:,1)),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro');hold on;
% h11=ploterr(log10(datalowall(:,1)),datalowall(:,3),{xbrd(:,1),xbrd(:,2)},datalowall(:,4),'ms','logx');hold on;
% h11=ploterr(log10(datacomp(:,1)),datacomp(:,3),[],datacomp(:,4),'gd');hold on;
% h11=ploterr(log10(datapassive(:,1)),datapassive(:,3),{xbrd(:,1),xbrd(:,2)},datapassive(:,4),'rd','logx');hold on;
% h11=ploterr(log10(dataZ02passive(:,1)),dataZ02passive(:,3),{xbrd(:,1),xbrd(:,2)},dataZ02passive(:,4),'rd','logx');hold on;
% h111=ploterr(log10(datapassive1R03offset(:,1)),datapassive1R03offset(:,3),[],[],'g-','logx');hold on;
% plot(log10(datapassive1R(:,1)),datapassive1R(:,3),'b-');hold on;
% h11=ploterr(log10(dataFlagSFRcomp(:,1)),dataFlagSFRcomp(:,3),{xbrd(:,1),xbrd(:,2)},dataFlagSFRcomp(:,4),'ro','logx');hold on;
h11=ploterr(log10(dataField(:,1)),dataField(:,3),{xbrd(:,1),xbrd(:,2)},dataField(:,4),'ro');hold on;
f=datahigh(:,1)<0.5e11; %only show Nlens>100 bins
% h111=ploterr(log10(datahighall(f,1)),datahighall(f,3),{xbrd(f,1),xbrd(f,2)},datahighall(f,4),'gd','logx');hold on;
% h111=ploterr(log10(dataactive(f,1)),dataactive(f,3),{xbrd(f,1),xbrd(f,2)},dataactive(f,4),'bs','logx');hold on;
% h111=ploterr(log10(dataactive1R(f,1)),dataactive1R(f,3),{xbrd(f,1),xbrd(f,2)},dataactive1R(f,4),'cx','logx');hold on;
% h111=ploterr(log10(dataactive1R03offset(f,1)),dataactive1R03offset(f,3),{xbrd(f,1),xbrd(f,2)},dataactive1R03offset(f,4),'gx','logx');hold on;
% h111=ploterr(log10(dataactiveNewEstimate(f,1)),dataactiveNewEstimate(f,3),{xbrd(f,1),xbrd(f,2)},dataactiveNewEstimate(f,4),'bs','logx');hold on;
h111=ploterr(log10(dataFieldHigh(f,1)),dataFieldHigh(f,3),{xbrd(f,1),xbrd(f,2)},dataFieldHigh(f,4),'gs');hold on;
% h=ploterr(log10(datahighz(:,1)),datahighz(:,3),[],[],'bx');set(h(1),'markersize',15);

% plot(x,15.2+1.58*(x-12),'b-');
% plot(x,15.01+1.43*(x-12),'r-');
% plot(x,13.96+0.93*(x-12),'g-');
% plot(x,13.75+0.78*(x-12),'g--');
% plot(x,13.59+0.85*(x-12),'r-');
xlim([9.5,12]);
ylim([10,15.5]);
% xscale('log');
xlabel('$\log(M_\star$[M$_\odot/h^2$])');
ylabel('$\log(M_h$[M$_\odot/h$])');
l=legend([h11(1),h111(1), h2,h22],'Field Central','Field high SSFR','Mock Central','Mock Active');
% l=legend([h1(1),h111(1),h11(1),h2,h22],'Combined','Active Central','Passive Central','Mock Central','Mock Active','$\sigma=0.2$ HOD');
% l=legend([h1(1),h2,h22,h222,h4,h5,h6,h7,h3(2)],'Data','Mock','Mock Active','Mock Passive','WangL13','WangLY13','Moster13','Guo10','All $\sigma=0.2$');
set(l,'location','northwest','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/M-Mstar-central-FieldHighSSFR.eps');
%% MMGs
dataMMG=[7.41951e+09   1.42993e+09  11.8       0.3974                  
1.45489e+10   2.87274e+09  12.02      0.3314                  
2.54499e+10   3.56541e+09  12.09      0.3327                  
4.0711e+10   5.75661e+09  12.43      0.2112                  
6.52116e+10   9.13312e+09  12.5       0.2765                  
1.05041e+11   1.49824e+10  13.32      0.131                   
1.71179e+11   2.26226e+10  13.29      0.2337                  
2.69282e+11   2.76495e+10  5.173      nan                   
4.90697e+11   6.3415e+10  13.01      1.061
]; %MMG, all
dataMMGgroup=[7.36608e+09   1.45245e+09  12.58      0.3982                  
1.4653e+10   2.8622e+09  12.45      0.3831                  
2.57412e+10   3.66671e+09  11.22      2.613                   
4.11456e+10   5.76043e+09  12.37      0.3671                  
6.62136e+10   9.72961e+09  12.52      0.3607                  
1.05292e+11   1.50783e+10  13.46      0.1403                  
1.70374e+11   2.2126e+10  13.43      0.243                   
2.6964e+11   3.153e+10  13.42      0.6137                  
5.09329e+11   4.19028e+10  13.31      1.143]; %MMG, grouped (EnvLevel>0)
dataMMGhigh=[7.30232e+09   1.43844e+09  12.24      0.3508                  
1.44093e+10   2.88207e+09  12.51      0.3114                  
2.52676e+10   3.50485e+09  12.63      0.4073                  
3.96241e+10   5.56163e+09  13.18      0.2953                  
6.32232e+10   8.74816e+09  5.425      814.7                   
1.00774e+11   1.37185e+10  13.78      0.7485                  
1.58685e+11   1.27925e+10  14.52      0.6313                  
2.61078e+11   2.69671e+10  10.28      25.58 
nan nan nan nan]; %MMG, highSSFR
dataMMGgrouphigh=[7.2291e+09   1.41912e+09  13.09      0.414                   
1.44677e+10   2.7834e+09  12.97      0.4839                  
2.57917e+10   3.6171e+09  11.91      2.306                   
4.03134e+10   5.56071e+09  4.176      nan                   
6.41118e+10   9.35455e+09  6.72       nan                   
1.00791e+11   1.39364e+10  14.4       0.4905                  
1.52695e+11   5.43435e+09  15.03      0.781 
nan nan nan nan]; %MMG, highSSFR, grouped
dataMMGgroupHighZ02=[7.16182e+09   1.41343e+09  13.35      0.3188                  
1.41703e+10   2.71338e+09  12.8       0.6719                  
2.53337e+10   3.46467e+09  12.73      0.9675                  
3.62463e+10   3.67427e+09  5.656      nan                   
6.3328e+10   7.92678e+09  6.141      nan                    
1.16574e+11   6.58097e+09  14.74      0.5914                  
1.58129e+11   0  15.48      0.8361
nan nan nan nan]; %MMG, highSSFR, grouped, z<0.2
dataGroupHigh=[7.21794e+09   1.4091e+09  11.51      2.296                   
1.45545e+10   2.80457e+09  12.92      0.484                   
2.56843e+10   3.64784e+09  12.99      0.5379                  
4.08173e+10   5.63286e+09  4.601      nan%907.4  216groups                  
6.40029e+10   9.36386e+09  5.822      nan %955.7                   
9.9929e+10   1.3537e+10  14.16      0.6667                  
1.4726e+11   0  8.424      nan%266.2 
nan nan nan nan]; %Group Centrals (EnvLevel>0 and CenSampleIter>0), high SSFR
dataGroup=[7.36856e+09   1.439e+09  10.81      5.043                   
1.49508e+10   2.94332e+09  12.38      0.445                   
2.59873e+10   3.69491e+09  11.78      1.282                   
4.21248e+10   5.78444e+09  12.38      0.3562                  
6.72419e+10   9.53521e+09  12.73      0.2806                  
1.07723e+11   1.50815e+10  13.38      0.1597                  
1.73102e+11   2.32395e+10  13.51      0.2278                  
2.74603e+11   3.73701e+10  13.65      0.4133                  
4.48152e+11   5.58734e+10  13.17      1.574];%All groups
dataGroupZ02=[7.33998e+09   1.43439e+09  11.74      1.261                   
1.4886e+10   2.95764e+09  12.29      0.5198                  
2.57734e+10   3.64964e+09  11.83      1.181                   
4.16951e+10   5.86541e+09  12.37      0.4069                  
6.63372e+10   9.3968e+09  12.9       0.2671                  
1.05858e+11   1.46249e+10  13.68      0.14                    
1.67697e+11   2.09796e+10  13.64      0.2942                  
2.72926e+11   3.23325e+10  14.35      0.2843                  
4.33198e+11   0  8.292      113.3    ]; %z<0.2 All groups (N>=2)
dataGroupHighZ02=[7.13093e+09   1.38362e+09  12.32      0.9711                  
1.41634e+10   2.72666e+09  12.91      0.5687                  
2.52212e+10   3.66758e+09  13.2       0.5971                  
3.68191e+10   3.82643e+09  4.593      nan %911.4                   
5.7485e+10   3.03102e+09  10.16      nan %28.48                   
1.22528e+11   0  14.99      0.6348       
nan nan nan nan
nan nan nan nan]; %z<0.2 Group Centrals (EnvLevel>0 and CenSampleIter>0), high SSFR
dataIterMMG=[7.42108e+09   1.43158e+09  12.33      0.6195                  
1.49894e+10   2.96746e+09  12.31      0.5538                  
2.60444e+10   3.68019e+09  6.355      nan                   
4.20871e+10   5.79628e+09  12.4       0.3754                  
6.75196e+10   9.62247e+09  12.59      0.3632                  
1.08012e+11   1.50716e+10  13.4       0.1627                  
1.7315e+11   2.30023e+10  13.44      0.2564                  
2.74865e+11   3.73038e+10  13.55      0.4725                  
4.50555e+11   5.58071e+10  13.05      1.748]; %Iter and MMG at the same time
dataIterMMGZ02=[7.40193e+09   1.43587e+09  12.5       0.5029                  
1.49286e+10   2.97615e+09  12.17      0.6643                  
2.58423e+10   3.64048e+09  11.73      1.43                    
4.16503e+10   5.81607e+09  12.37      0.4355                  
6.63876e+10   9.44048e+09  12.9       0.2784                  
1.06105e+11   1.45552e+10  13.71      0.1391                  
1.68073e+11   2.1283e+10  13.57      0.3231                  
2.72926e+11   3.23325e+10  14.35      0.2843                  
4.33198e+11   0  8.292      nan]; %z<0.2 IterMMG
dataIterMMGHighZ02=[7.15077e+09   1.40635e+09  13.05      0.495                   
1.43049e+10   2.74403e+09  11.21      6.099                   
2.54074e+10   3.56948e+09  12.92      0.9042                  
3.64149e+10   3.68488e+09  4.347      nan                    
5.7485e+10   3.03102e+09  10.16      nan                   
1.22528e+11   0  14.99      0.6348
nan nan nan nan
nan nan nan nan
nan nan nan nan]; %z<0.2 IterMMG and high SSFR
dataIterMMGHigh=[7.23833e+09   1.4098e+09  12.7       0.7028                  
1.45395e+10   2.77048e+09  12.48      1.322                   
2.58801e+10   3.61487e+09  12.45      1.334                   
4.06796e+10   5.61686e+09  4.289      nan                    
6.42721e+10   9.77489e+09  7.906      nan                   
1.00451e+11   1.4685e+10  14.58      0.4777                  
1.4726e+11   0  8.424      266.2 
nan nan nan nan
nan nan nan nan]; %IterMMG and high SSFR
dataIterNoneMMGZ02=[7.2128e+09   1.42293e+09  5.021      nan                   
1.47067e+10   2.87152e+09  12.7       0.7322                  
2.51883e+10   3.67478e+09  12.46      1.409                   
4.21469e+10   6.32464e+09  12.39      1.15                    
6.53949e+10   8.4831e+09  12.92      0.9373                  
1.0045e+11   1.50939e+10  13         1.393                   
1.58567e+11   6.51303e+09  14.77      0.4451
nan nan nan nan
nan nan nan nan];%IterNoneMMG, z<0.2
dataIterNoneMMGHighZ02=[7.2128e+09   1.42293e+09  5.048      nan                   
1.47067e+10   2.87152e+09  13.84      0.3563                  
2.51883e+10   3.67478e+09  13.72      0.6264                  
4.21469e+10   6.32464e+09  14.15      1.068                   
6.53949e+10   8.4831e+09  nan nan                       
1.0045e+11   1.50939e+10  nan         nan                       
1.58567e+11   6.51303e+09  nan         nan  
nan nan nan nan; nan nan nan nan;];%IterNoneMMG, highSSFR, z<0.2


xbrd=[5e9,1e10; 1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11; 3.7e11, 6e11];
xbrd=log10(xbrd);
myfigure;

%-------standard mock--------------
f=galmock.CentralSampleSMpeak>0&galmock.EnvLevel>0&galmock.Z<0.2; %about 90% are real centrals
x=linspace(8,12,12);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(galmock.Mstar(f)),log10(galmock.Mhalo(f)),x,0.68);
h2=plot(xmed,ymed,'-','linewidth',2,'color','r');hold on;
% plot(xmed,yl,'r.-');
f=f&galmock.SFR./galmock.Mstar>10^-1.5;
[xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(galmock.Mstar(f)),log10(galmock.Mhalo(f)),x,0.68);
h22=plot(xmed,ymed,'-','linewidth',2,'color','b');hold on;

% h0=plot(log10(dataIterNoneMMGZ02(:,1)),dataIterNoneMMGZ02(:,3),'m--');
% h0=ploterr(log10(dataIterNoneMMGZ02(:,1)),dataIterNoneMMGZ02(:,3),[],dataIterNoneMMGZ02(:,4),'ms--');hold on;
h11=ploterr(log10(dataMMGgroup(:,1)),dataMMGgroup(:,3),{xbrd(:,1),xbrd(:,2)},dataMMGgroup(:,4),'ro');hold on;
% f=dataMMGhigh(:,1)<0.5e11; %only show Nlens>100 bins
% h111=ploterr(log10(dataMMGhigh(f,1)),dataMMGhigh(f,3),{xbrd(f,1),xbrd(f,2)},dataMMGhigh(f,4),'gs');hold on;
% h11=ploterr(log10(dataMMGgroup(:,1)),dataMMGgroup(:,3),{xbrd(:,1),xbrd(:,2)},dataMMGgroup(:,4),'ro');hold on;
f=dataMMGhigh(:,1)<0.5e13; %only show Nlens>100 bins
% h111=ploterr(log10(dataMMGgrouphigh(f,1)),dataMMGgrouphigh(f,3),{xbrd(f,1),xbrd(f,2)},dataMMGgrouphigh(f,4),'gs');hold on;
% h11=ploterr(log10(dataGroup(:,1)),dataGroup(:,3),{xbrd(:,1),xbrd(:,2)},dataGroup(:,4),'mo');hold on;
f=dataGroupHighZ02(:,1)<inf; %only show Nlens>100 bins
% h11=ploterr(log10(dataGroupHighZ02(f,1)),dataGroupHighZ02(f,3),{xbrd(f,1),xbrd(f,2)},dataGroupHighZ02(f,4),'bs');hold on;
h111=ploterr(log10(dataMMGgroupHighZ02(f,1)),dataMMGgroupHighZ02(f,3),{xbrd(f,1),xbrd(f,2)},dataMMGgroupHighZ02(f,4),'go');hold on;
% h11=ploterr(log10(dataIterNoneMMGHighZ02(f,1)),dataIterNoneMMGHighZ02(f,3),{xbrd(f,1),xbrd(f,2)},dataIterNoneMMGHighZ02(f,4),'cx');hold on;
h0=plot(log10(dataIterMMGZ02(f,1)),dataIterMMGZ02(f,3),'m--');hold on;
% h0=ploterr(log10(dataGroupZ02(f,1)),dataGroupZ02(f,3),{xbrd(f,1),xbrd(f,2)},dataGroupZ02(f,4),'kd');hold on;

xlim([9.5,12]);
ylim([10,15.5]);
% xscale('log');
xlabel('$\log(M_\star$[M$_\odot/h^2$])');
ylabel('$\log(M_h$[M$_\odot/h$])');
l=legend([h11(1),h111(1),h0(1),h2,h22],'GroupMMG', 'GroupMMG active','IterMMG','MockGroupMMG','MockGroupMMG Active');
% l=legend([h1(1),h111(1),h11(1),h2,h22],'Combined','Active Central','Passive Central','Mock Central','Mock Active','$\sigma=0.2$ HOD');
% l=legend([h1(1),h2,h22,h222,h4,h5,h6,h7,h3(2)],'Data','Mock','Mock Active','Mock Passive','WangL13','WangLY13','Moster13','Guo10','All $\sigma=0.2$');
set(l,'location','southeast','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/SatContamination/M-Mstar-MMGgroup-VolLimited.eps');
%% M=A(M*) comparison, z dep             
data=[7.3685e+09   1.45222e+09  11.51      0.5988                  
1.4683e+10   2.90974e+09  11.97      0.3629                  
2.57841e+10   3.59833e+09  12.16      0.3087                  
4.16434e+10   5.83513e+09  12.43      0.2122                  
6.67911e+10   9.34632e+09  12.61      0.2369                  
1.07375e+11   1.52612e+10  13.26      0.1444                  
1.73309e+11   2.39779e+10  13.34      0.2277                  
2.77191e+11   3.80254e+10  nan      nan                   
4.48464e+11   5.69671e+10  12.9       1.273
];%r<19.4 centrals
xbrd=[5e9,1e10; 1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11; 3.7e11, 6e11];
xbrd=log10(xbrd);
datahigh=[7.29997e+09   1.43735e+09  12.01      0.4588                  
1.442e+10   2.88151e+09  12.51      0.3106                  
2.52676e+10   3.51158e+09  12.8       0.3169                  
3.97499e+10   5.5919e+09  13.11      0.3219                  
6.31942e+10   8.748e+09  5.297     nan                   %266 halos
1.00395e+11   1.35373e+10  13         nan                   %41 halos
1.5887e+11   1.47669e+10  14.22      1.096                   %3 halos in total
2.61078e+11   2.69671e+10  10.28      nan                    %3 halos in total
nan nan nan nan
]; %high SSFR, and restricted to r<19.4
datahighz=[7.30049e+09   1.45166e+09  11.96      0.6161                  
1.46217e+10   2.87625e+09  12.49      0.3783                  
2.59301e+10   3.58184e+09  12.23      0.7441                  
4.15344e+10   5.82956e+09  12.64      0.3616                  
6.67339e+10   9.31684e+09  12.69      0.572                   
1.07675e+11   1.54596e+10  12.46      1.126                   
1.74241e+11   2.41788e+10  13.53      0.4475                  
2.76378e+11   3.8721e+10  8.453      nan                  
4.39403e+11   5.11245e+10  13.25      1.325 ];%high-z subsample
myfigure;
sigma_common=0.2;
for i=1:numel(ms)
    zwang(i)=medianMh(log10(ms(i)),sigma_common,Pwang4);
    zling(i)=medianMh(log10(ms(i)),sigma_common,Pling);
    zmosterZ0(i)=medianMh(log10(ms(i)),sigma_common,PmosterZ0);
    zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo); 
end
zbound=minmax([zwang',zling',zmosterZ0',zguo']);
h3=area(log10(ms),[zbound(:,1),zbound(:,2)-zbound(:,1)]);hold on;
set(h3(2),'facecolor','y','edgecolor','w');
set(h3(1),'facecolor','w','edgecolor','w');
set(gca,'layer','top');
%----------mock predictions----------------
x=linspace(8,12,12);
%-------standard mock--------------
f=galmock.CentralSampleIter>0; %about 90% are real centrals
[xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(galmock.Mstar(f)),log10(galmock.Mhalo(f)),x,0.68);
h2=plot(xmed,ymed,'-','linewidth',2,'color','r');hold on;
% end
% plot(xmed,yl,'r.-');
f=galmock.CentralSampleIter>0&galmock.Z>0.0228*log10(galmock.Mstar).^2-0.341*log10(galmock.Mstar)+1.307;
[xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(galmock.Mstar(f)),log10(galmock.Mhalo(f)),x,0.68);
h22=plot(xmed,ymed,'-','linewidth',2,'color','b');hold on;
% plot(xmed,yl,'--');

h1=ploterr(log10(data(:,1)),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro');hold on;
h11=ploterr(log10(datahighz(:,1)),datahighz(:,3),{xbrd(:,1),xbrd(:,2)},datahighz(:,4),'gd');hold on;
xlim([9.5,12]);
ylim([10,15.5]);
% xscale('log');
xlabel('$\log(M_\star$[M$_\odot/h^2$])');
ylabel('$\log(M_h$[M$_\odot/h$])');
l=legend([h1(1),h11(1),h2,h22,h3(2)],'All Central','High-z Central','Mock Central','Mock High-z','$\sigma=0.2$ HOD');
% l=legend([h1(1),h2,h22,h222,h4,h5,h6,h7,h3(2)],'Data','Mock','Mock Active','Mock Passive','WangL13','WangLY13','Moster13','Guo10','All $\sigma=0.2$');
set(l,'location','northwest','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/M-Mstar-central-Zhigh.eps');
%% zbins            
data=[7.3685e+09   1.45222e+09  11.51      0.5988                  
1.4683e+10   2.90974e+09  11.97      0.3629                  
2.57841e+10   3.59833e+09  12.16      0.3087                  
4.16434e+10   5.83513e+09  12.43      0.2122                  
6.67911e+10   9.34632e+09  12.61      0.2369                  
1.07375e+11   1.52612e+10  13.26      0.1444                  
1.73309e+11   2.39779e+10  13.34      0.2277                  
2.77191e+11   3.80254e+10  nan      nan                   
4.48464e+11   5.69671e+10  12.9       1.273
];%r<19.4 centrals
xbrd=[5e9,1e10; 1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11; 3.7e11, 6e11];
xbrd=log10(xbrd);
data1=[7.18464e+09   1.4273e+09  10.52      4.331                   
1.41736e+10   2.82257e+09  12.27      0.5071                  
2.50457e+10   3.64081e+09  12.44      0.5452                  
4.14283e+10   5.9643e+09  12.15      1.09                    
6.60571e+10   9.86382e+09  13.09      0.3901                  
1.05064e+11   1.46322e+10  13.49      0.4159                  
1.95752e+11   1.41766e+10  14.74      0.3863                  
2.8812e+11   3.24175e+10  13.57      1.283                   
4.33198e+11   0  8.292      nan    
]; %z=0-0.1
data2=[7.30834e+09   1.4559e+09  11.86      0.4645                  
1.44243e+10   2.90575e+09  11.3       1.337                   
2.53026e+10   3.54066e+09  12.06      0.4739                  
4.10152e+10   5.74464e+09  12.53      0.2748                  
6.54969e+10   9.07545e+09  12.7       0.3118                  
1.05008e+11   1.43576e+10  13.67      0.1387                  
1.66913e+11   2.0985e+10  13.51      0.3336                  
2.72216e+11   2.86494e+10  6.503      nan                   
5.86511e+11   0  10.37      nan      ]; %z=0.1-0.2
data3=[7.65118e+09   1.41472e+09  5.956      nan                   
1.50754e+10   2.86852e+09  12.55      0.3622                  
2.60748e+10   3.59041e+09  12.27      0.4867                  
4.1536e+10   5.77561e+09  12.2       0.4588                  
6.5994e+10   9.19941e+09  12.49      0.473                   
1.05349e+11   1.48097e+10  12.62      0.6434                  
1.70447e+11   2.29394e+10  13.53      0.278                   
2.73411e+11   3.42799e+10  12.34      4.568                   
4.60995e+11   6.44675e+10  13.71      0.8511];%z=0.2-0.3
datahighz=[7.30049e+09   1.45166e+09  11.96      0.6161                  
1.46217e+10   2.87625e+09  12.49      0.3783                  
2.59301e+10   3.58184e+09  12.23      0.7441                  
4.15344e+10   5.82956e+09  12.64      0.3616                  
6.67339e+10   9.31684e+09  12.69      0.572                   
1.07675e+11   1.54596e+10  12.46      1.126                   
1.74241e+11   2.41788e+10  13.53      0.4475                  
2.76378e+11   3.8721e+10  8.453      nan                  
4.39403e+11   5.11245e+10  13.25      1.325 ];%high-z subsample
myfigure;
sigma_common=0.2;
for i=1:numel(ms)
    zwang(i)=medianMh(log10(ms(i)),sigma_common,Pwang4);
    zling(i)=medianMh(log10(ms(i)),sigma_common,Pling);
    zmosterZ0(i)=medianMh(log10(ms(i)),sigma_common,PmosterZ0);
    zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo); 
end
zbound=minmax([zwang',zling',zmosterZ0',zguo']);
h3=area(log10(ms),[zbound(:,1),zbound(:,2)-zbound(:,1)]);hold on;
set(h3(2),'facecolor','y','edgecolor','w');
set(h3(1),'facecolor','w','edgecolor','w');
set(gca,'layer','top');
%----------mock predictions----------------
x=linspace(8,12,12);
%-------standard mock--------------
f=galmock.CentralSampleIter>0&galmock.Z<0.1; %about 90% are real centrals
[xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(galmock.Mstar(f)),log10(galmock.Mhalo(f)),x,0.68);
h2=plot(xmed,ymed,'-','linewidth',2,'color','r');hold on;
f=galmock.CentralSampleIter>0&galmock.Z>0.1&galmock.Z<0.2; %about 90% are real centrals
[xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(galmock.Mstar(f)),log10(galmock.Mhalo(f)),x,0.68);
h22=plot(xmed,ymed,'-','linewidth',2,'color','g');hold on;
f=galmock.CentralSampleIter>0&galmock.Z>0.2&galmock.Z<0.3; %about 90% are real centrals
[xmed,ymed,yl,xm,ym,yerr]=skeleton(log10(galmock.Mstar(f)),log10(galmock.Mhalo(f)),x,0.68);
h222=plot(xmed,ymed,'-','linewidth',2,'color','b');hold on;
% h1=ploterr(log10(data(:,1)),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro');hold on;
h1=ploterr(log10(data1(:,1)),data1(:,3),{xbrd(:,1),xbrd(:,2)},data1(:,4),'ro');hold on;
h11=ploterr(log10(data2(:,1)),data2(:,3),{xbrd(:,1),xbrd(:,2)},data2(:,4),'gd');hold on;
h111=ploterr(log10(data3(:,1)),data3(:,3),{xbrd(:,1),xbrd(:,2)},data3(:,4),'bs');hold on;
% h11=ploterr(log10(datahighz(:,1)),datahighz(:,3),{xbrd(:,1),xbrd(:,2)},datahighz(:,4),'gd');hold on;
xlim([9.5,12]);
ylim([10,15.5]);
% xscale('log');
xlabel('$\log(M_\star$[M$_\odot/h^2$])');
ylabel('$\log(M_h$[M$_\odot/h$])');
l=legend([h1(1),h11(1),h111(1),h2,h22,h222,h3(2)],'z=0-0.1','0.1-0.2','0.2-0.3','Mock 0-0.1','Mock 0.1-0.2','Mock 0.2-0.3','$\sigma=0.2$ HOD');
% l=legend([h1(1),h2,h22,h222,h4,h5,h6,h7,h3(2)],'Data','Mock','Mock Active','Mock Passive','WangL13','WangLY13','Moster13','Guo10','All $\sigma=0.2$');
set(l,'location','northwest','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/M-Mstar-central-Zbins.eps');
%% color dep             
data=[7.3685e+09   1.45222e+09  11.51      0.5988                  
1.4683e+10   2.90974e+09  11.97      0.3629                  
2.57841e+10   3.59833e+09  12.16      0.3087                  
4.16434e+10   5.83513e+09  12.43      0.2122                  
6.67911e+10   9.34632e+09  12.61      0.2369                  
1.07375e+11   1.52612e+10  13.26      0.1444                  
1.73309e+11   2.39779e+10  13.34      0.2277                  
2.77191e+11   3.80254e+10  nan      nan                   
4.48464e+11   5.69671e+10  12.9       1.273
];%r<19.4 centrals
xbrd=[5e9,1e10; 1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11; 3.7e11, 6e11];
datablue=[7.26127e+09   1.43877e+09  11.31      0.9383                  
1.43845e+10   2.89164e+09  11.99      0.5554                  
2.54431e+10   3.56093e+09  11.81      1.145                   
4.08797e+10   5.77946e+09  12.5       0.4886                  
6.62334e+10   9.37552e+09  4.911      nan                   
1.08646e+11   1.54148e+10  6.302      nan                   
1.74113e+11   2.49938e+10  13.64      0.4852                  
2.80732e+11   3.71322e+10  7.116      nan                   
4.6118e+11   6.19565e+10  11.25      nan]; %blue
datared=[7.6439e+09   1.45034e+09  11.81      0.7046                  
1.50055e+10   2.89478e+09  11.95      0.4776                  
2.59673e+10   3.60498e+09  12.24      0.3086                  
4.18916e+10   5.83162e+09  12.41      0.2384                  
6.69035e+10   9.33639e+09  12.67      0.2249                  
1.07125e+11   1.52182e+10  13.35      0.1336                  
1.73136e+11   2.37504e+10  13.3       0.2519                  
2.76161e+11   3.82198e+10  8.034      nan                   
4.44404e+11   5.46594e+10  13.3       1.022 ]; %red
databluehigh=[7.2836e+09   1.43616e+09  11.91      0.5504                  
1.42952e+10   2.86541e+09  12.59      0.3139                  
2.49872e+10   3.41752e+09  12.8       0.3918                  
3.94474e+10   5.47259e+09  13.27      0.3479                  
6.31573e+10   8.71878e+09  13.44      0.7397                  
9.52068e+10   9.99711e+09  8.435      nan                   
1.64675e+11   1.50332e+10  14.27      1.03                    
2.78027e+11   1.51316e+10  10.34      nan
nan nan nan nan]; %blue and high
myfigure;
sigma_common=0.2;
for i=1:numel(ms)
    zwang(i)=medianMh(log10(ms(i)),sigma_common,Pwang4);
    zling(i)=medianMh(log10(ms(i)),sigma_common,Pling);
    zmosterZ0(i)=medianMh(log10(ms(i)),sigma_common,PmosterZ0);
    zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo); 
end
zbound=minmax([zwang',zling',zmosterZ0',zguo']);
h3=area(ms,[zbound(:,1),zbound(:,2)-zbound(:,1)]);hold on;
set(h3(2),'facecolor','y','edgecolor','w');
set(h3(1),'facecolor','w','edgecolor','w');
set(gca,'layer','top');
%----------mock predictions----------------
x=logspace(8,12,12);
%-------standard mock--------------
f=galmock.CentralSampleIter>0; %about 90% are real centrals
[xmed,ymed,yl,xm,ym,yerr]=skeleton(galmock.Mstar(f),log10(galmock.Mhalo(f)),x,0.68);
% h2=plot(xmed,ymed,'-','linewidth',2,'color','r');hold on;
% end
% plot(xmed,yl,'r--');
f=galmock.CentralSampleIter>0&galmock.color<1.6;
[xmed,ymed,yl,xm,ym,yerr]=skeleton(galmock.Mstar(f),log10(galmock.Mhalo(f)),x,0.68);
h22=plot(xmed,ymed,'-','linewidth',2,'color','b');hold on;
plot(xmed,yl,'b--');
f=galmock.CentralSampleIter>0&galmock.color>1.6;
x=logspace(8,12,15);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(galmock.Mstar(f),log10(galmock.Mhalo(f)),x,0.68);
h222=plot(xmed,ymed,'-','linewidth',2,'color','r');hold on;
plot(xmed,yl,'r--');
% ------------volume limited mock----------------
% magcut=19.4; %almost no effect in applying a flux limit or not, in the mocks.
% tmpmock=guo;
% f=tmpmock.rmag<magcut;
% [xmed,ymed]=skeleton(tmpmock.sm(f),log10(tmpmock.mh(f,1)),x,0.68);
% % h220=plot(xmed,ymed,'--','linewidth',2,'color','r');hold on;
% f=tmpmock.rmag<magcut&tmpmock.color<2;
% [xmed,ymed]=skeleton(tmpmock.sm(f),log10(tmpmock.mh(f,1)),x,0.68);
% h22=plot(xmed,ymed,'--','linewidth',2,'color','b');hold on;
% f=tmpmock.rmag<magcut&tmpmock.color>2;
% [xmed,ymed]=skeleton(tmpmock.sm(f),log10(tmpmock.mh(f,1)),x,0.68);
% h222=plot(xmed,ymed,'--','linewidth',2,'color','m');hold on;
% tmpmock=guo;
% f=tmpmock.rmag<magcut;
% [xmed,ymed]=skeleton(tmpmock.sm(f),log10(tmpmock.mh(f,1)),x,0.68);
% h220=plot(xmed,ymed,'p--','linewidth',2,'color','r');hold on;
% f=tmpmock.rmag<magcut&tmpmock.sfr./tmpmock.sm>10^-1.5&tmpmock.z<0.31;
% [xmed,ymed]=skeleton(tmpmock.sm(f),log10(tmpmock.mh(f,1)),x,0.68);
% h220=plot(xmed,ymed,'p--','linewidth',2,'color','b');hold on;
% tmpmock=font;
% f=tmpmock.rmag<magcut;
% [xmed,ymed]=skeleton(tmpmock.sm(f),log10(tmpmock.mh(f,1)),x,0.68);
% h220=plot(xmed,ymed,'<--','linewidth',2,'color','r');hold on;
% f=tmpmock.rmag<magcut&tmpmock.sfr./tmpmock.sm>10^-1.5&tmpmock.z<0.31;
% [xmed,ymed]=skeleton(tmpmock.sm(f),log10(tmpmock.mh(f,1)),x,0.68);
% h220=plot(xmed,ymed,'<--','linewidth',2,'color','b');hold on;

% h1=ploterr(data(:,1),data(:,3),{xbrd(:,1),xbrd(:,2)},data(:,4),'ro','logx');hold on;
h11=ploterr(datared(:,1),datared(:,3),{xbrd(:,1),xbrd(:,2)},datared(:,4),'md','logx');hold on;
% h11=ploterr(datared(:,1),datared(:,3),[],datared(:,4),'md','logx');hold on;
% h11=ploterr(datazbound(:,1),datazbound(:,3),{xbrd(:,1),xbrd(:,2)},datazbound(:,4),'gx','logx');hold on;
f=datablue(:,1)<1e12;
h111=ploterr(datablue(f,1),datablue(f,3),{xbrd(f,1),xbrd(f,2)},datablue(f,4),'bs','logx');hold on;
ploterr(databluehigh(:,1),databluehigh(:,3),[],databluehigh(:,4),'gs','logx');hold on;
xlim([4e9,1e12]);
ylim([10,15.5]);
xscale('log');
xlabel('$M_\star$[M$_\odot/h^2$]');
ylabel('$\log(M_h$[M$_\odot/h$])');
l=legend([h111(1),h11(1),h22,h222,h3(2)],'Blue Central','Red Central','Mock Blue','Mock Red','$\sigma=0.2$ HOD');
% l=legend([h1(1),h2,h22,h222,h4,h5,h6,h7,h3(2)],'Data','Mock','Mock Active','Mock Passive','WangL13','WangLY13','Moster13','Guo10','All $\sigma=0.2$');
set(l,'location','northwest','interpreter','latex');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/M-Mstar-central-color-guo.eps');
%% sample distribution
nbin=30;
smth=@(n) conv2(n,[0.05,0.1,0.05;0.1,0.4,0.1;0.05,0.1,0.05;],'same');
% smth=@(n) conv2(n,[0.08,0.15,0.08;0.15,0.23,0.15;0.08,0.15,0.08;],'same');
percents=[0.3,0.6,0.9];
myfigure;
f=gal.FlagSFR>0;
x=log10(gal.SMsps(f));y=log10(gal.SFR(f)./gal.SMsps(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[-5,4]);
[~,h0]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'r');hold on;%colormap(cmap);
f=galmock.zhub<0.31;
x=log10(galmock.Mstar(f));y=log10(galmock.SFR(f)./galmock.Mstar(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[-5,4]);
[~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'g');hold on;%colormap(cmap);
% f=rawmock.z<0.31;
% x=log10(rawmock.sm(f));y=log10(rawmock.SFR(f)./rawmock.sm(f));
% [xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[-5,4]);
% [~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'b');hold on;%colormap(cmap);
plot([8,12],[-1.5,-1.5],'k--');
xlabel('log(Mstar)');ylabel('log(SSFR)');
l=legend('Data','Mock');set(l,'interpreter','latex');
xlabel('$\log(M_\star$[M$_\odot/h^2$])');
ylabel('$\log$(SSFR[Gyr$^{\frac{\ }{\ }1}$])');
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/SSFR-Mstar.eps');
%% SSFR-M*/L
RmagSun=4.67;
UmagSun=6.77;
nbin=30;
smth=@(n) conv2(n,[0.05,0.1,0.05;0.1,0.4,0.1;0.05,0.1,0.05;],'same');
% smth=@(n) conv2(n,[0.08,0.15,0.08;0.15,0.23,0.15;0.08,0.15,0.08;],'same');
percents=[0.3,0.6,0.9];
myfigure;
f=gal.FlagSFR>0;
x=log10(gal.SMsps(f))+0.4*(gal.Umodel_abs(f)-RmagSun);y=log10(gal.SFR(f)./gal.SMsps(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[-2,2],[-5,4]);
[~,h0]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'r');hold on;%colormap(cmap);
f=galmock.zhub<0.31;
x=log10(galmock.Mstar(f))+0.4*(galmock.Umodel_abs(f)-RmagSun);y=log10(galmock.SFR(f)./galmock.Mstar(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[-2,2],[-5,4]);
[~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'g');hold on;%colormap(cmap);
% f=rawmock.z<0.31;
% x=log10(rawmock.sm(f));y=log10(rawmock.SFR(f)./rawmock.sm(f));
% [xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[-5,4]);
% [~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'b');hold on;%colormap(cmap);
plot([-2,2],[-1.5,-1.5],'k--');
xlabel('log(Mstar)');ylabel('log(SSFR)');
legend('Data','Mock');
xlabel('$\log(M_\star$[M$_\odot/h^2$])');
ylabel('$\log$(SSFR[Gyr$^{\frac{\ }{\ }1}$])');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/SSFR-Mstar.eps');
%%
nbin=30;
smth=@(n) conv2(n,[0.05,0.1,0.05;0.1,0.4,0.1;0.05,0.1,0.05;],'same');
% smth=@(n) conv2(n,[0.08,0.15,0.08;0.15,0.23,0.15;0.08,0.15,0.08;],'same');
percents=0.15:0.2:0.95;
f=gal.FlagSFR;
x=log10(gal.SMsps(f));y=log10(gal.SFR(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[-5,15]);
[~,h0]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'r');hold on;%colormap('jet');
f=galmock.zhub<0.31;
x=log10(galmock.Mstar(f));y=log10(galmock.SFR(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[-5,15]);
[~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'g');hold on;%colormap('jet');
xlabel('log(Mstar)');ylabel('log(SFR)');
legend('Data','Mock');
%%
f=galmock.zhub<0.31;
x=log10(galmock.Mstar(f));y=log10(galmock.SFR(f));
plot(x,y,'g.','markersize',1); hold on;
f=gal.FlagSFR;
x=log10(gal.SMsps(f));y=log10(gal.SFR(f));
plot(x,y,'r.','markersize',1);
hold on;
%%
myfigure;
y=log10(gal.SFR);
f=gal.IsAGN==0&gal.CentralSampleIter>0&y<gal.SMsps<10^10.5;
linhist(y(f),7.5:0.1:11,'stairsnorm','r-');hold on;
f=galmock.CentralSampleIter>0&galmock.Z<0.31&galmock.Mstar<10.^10.5;
y=log10(galmock.SFR);
y=y+normrnd(0,0.37,size(y));
y=y+0.12;
linhist(y(f),7.5:0.1:11,'stairsnorm','g--');
%%
myfigure;
y=log10(gal.SMsps);
f=gal.IsAGN==0&gal.CentralSampleIter>0;
linhist(y(f),8:0.1:10.5,'stairsnorm','r-');hold on;
f=galmock.CentralSampleIter>0&galmock.Z<0.31;
y=log10(galmock.Mstar(f));
y=y+normrnd(0,0.43,size(y));
y=y+0.05;
linhist(y,7:0.1:10.5,'stairsnorm','g--');
%% Sample definitions
myfigure;
y=log10(gal.SSFR);
f=gal.FlagSFR>=2&gal.CentralSampleIter>0;
linhist(y(f),-3:0.1:1,'stairsnorm','r-');hold on;
f=gal.IsAGN==0&gal.CentralSampleIter>0;
% f=gal.FlagSFR>=2&gal.FlagBalmer==1&gal.CentralSampleIter>0;
linhist(y(f),-3:0.1:1,'stairsnorm','g-');hold on;
f=gal.FlagSFR>=3&gal.CentralSampleIter>0; %2/3 of all the SFcentrals
linhist(y(f),-3:0.1:1,'stairsnorm','b-');hold on;
f=gal.FlagSFR>=2&gal.CentralSampleIter==0; %sat sample
linhist(y(f),-3:0.1:1,'stairsnorm','k-');hold on;
f=galmock.CentralSampleIter>0&galmock.Z<0.31&galmock.Mstar<10.^10.5;
y=log10(galmock.SSFR(f));
[~,~,~,h]=linhist(y,-3:0.1:1,'stairsnorm','m--');
y=y+normrnd(0,0.4,size(y))+0.08;
[~,~,~,h]=linhist(y,-3:0.1:1,'stairsnorm','c--');
l=legend('Clean Central','SF Central','SF Clean Central','Clean Sat','Mock Central','Mock Noise-added');set(l,'interpreter','latex','location','northwest');
plot([-1.5,-1.5],[0,0.2],'k:');
xlabel('$\log$(SSFR[Gyr$^{\frac{\ }{\ }1}$])');
ylabel('Fraction');
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/HistSSFR-samples.eps')
%%
myfigure;
f=gal.FlagSFR;
f=f&(gal.EnvLevel==0|gal.IsIterCen);
x=log10(gal.SMsps(f));y=log10(gal.SFR(f)./gal.SMsps(f));
linhist(y,-3:0.1:1,'stairsnorm');hold on;
f=galmock.Z<0.31;
f=f&(galmock.EnvLevel==0|galmock.IsIterCen);
x=log10(galmock.Mstar(f));y=log10(galmock.SFR(f)./galmock.Mstar(f));
[~,~,~,h]=linhist(y,-3:0.1:1,'stairsnorm','r--');
l=legend('Data','Mock');set(l,'interpreter','latex');
plot([-1.5,-1.5],[0,0.2],'k:');
xlabel('$\log$(SSFR[Gyr$^{\frac{\ }{\ }1}$])');
ylabel('Fraction');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/HistSSFR.eps')
%% Rmag vs. Mstar
nbin=40;
smth=@(n) conv2(n,[0.05,0.1,0.05;0.1,0.4,0.1;0.05,0.1,0.05;],'same');
% smth=@(n) conv2(n,[0.08,0.15,0.08;0.15,0.23,0.15;0.08,0.15,0.08;],'same');
percents=[0.1,0.3,0.6,0.9];
myfigure;
f=gal.FlagSFR>0;
x=(gal.Rpetro(f));y=log10(gal.SFR(f)./gal.SMsps(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[15,20],[-5,4]);
[~,h0]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'r');hold on;%colormap(cmap);
f=galmock.zhub<0.31;
x=(galmock.r_mag(f));y=log10(galmock.SFR(f)./galmock.Mstar(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[15,20],[-5,4]);
[~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'g');hold on;%colormap(cmap);
% f=rawmock.z<0.31;
% x=log10(rawmock.sm(f));y=log10(rawmock.SFR(f)./rawmock.sm(f));
% [xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[-5,4]);
% [~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'b');hold on;%colormap(cmap);
plot([8,12],[-1.5,-1.5],'k--');
xlabel('log(Mstar)');ylabel('log(SSFR)');
legend('Data','Mock');
xlabel('$\log(M_\star$[M$_\odot/h^2$])');
ylabel('$\log$(SSFR[Gyr$^{\frac{\ }{\ }1}$])');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/SSFR-Mstar.eps');
%% SFR vs. Mstar
nbin=40;
smth=@(n) conv2(n,[0.05,0.1,0.05;0.1,0.4,0.1;0.05,0.1,0.05;],'same');
% smth=@(n) conv2(n,[0.08,0.15,0.08;0.15,0.23,0.15;0.08,0.15,0.08;],'same');
percents=[0.05,0.15,0.3,0.6,0.9, 0.95];
myfigure;
f=gal.FlagSFR>0;
x=log10(gal.SMsps(f));y=log10(gal.SFR(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[5,12]);
[~,h0]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'r');hold on;%colormap(cmap);
f=galmock.zhub<0.31;
x=log10(galmock.Mstar(f));y=log10(galmock.SFR(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[5,12]);
[~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'g');hold on;%colormap(cmap);
% f=rawmock.z<0.31;
% x=log10(rawmock.sm(f));y=log10(rawmock.SFR(f)./rawmock.sm(f));
% [xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[-5,4]);
% [~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'b');hold on;%colormap(cmap);
plot([5,12],[5,12]+[-1.5,-1.5],'k--');
xlabel('log(Mstar)');ylabel('log(SFR)');
legend('Data','Mock');
xlabel('$\log(M_\star$[M$_\odot/h^2$])');
ylabel('$\log$(SFR[M$_\odot/h^2$Gyr$^{\frac{\ }{\ }1}$])');
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/SFR-Mstar.eps');
%% stellar mass populations
myfigure;
[xm1,ym1]=loghist(gal.SMsps(gal.CentralSampleIter>0&(gal.FlagSFR>0&gal.SSFR>10^-1.5)),[xbrd(:,1);xbrd(end)])
[xm,ym]=loghist(gal.SMsps(gal.CentralSampleIter>0),[xbrd(:,1);xbrd(end)]);
plot(xm,ym1./ym,'r');hold on;
[xm1,ym1]=loghist(galmock.Mstar(galmock.CentralSampleIter>0&galmock.SSFR>10^-1.5&galmock.zhub<0.31),[xbrd(:,1);xbrd(end)]);
[xm,ym]=loghist(galmock.Mstar(galmock.CentralSampleIter>0),[xbrd(:,1);xbrd(end)]);
plot(xm,ym1./ym,'g');
xscale('log');
yscale('log');
xlabel('$M_\star[\rm{M}_\odot/h^2]$');
ylabel('Active Fraction');
legend('Data','Mock');
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/ActiveFraction.eps');
%% color distr
myfigure;
x=-1:0.05:4;
gal.color=gal.Umodel_abs-gal.Rmodel_abs;
linhist(gal.color,x,'norm','r')
hold on;
galmock.color=galmock.Umodel_abs-galmock.Rmodel_abs;
linhist(galmock.color,x,'norm','b')
legend('Data','Mock');
plot([2,2],[0,0.1],'r--');
plot([1.6,1.6],[0,0.1],'b--');
ylim([0,0.06]);
xlabel('U-R');
ylabel('Probability');
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/ColorDistr.eps');
%% SSFR-color
nbin=40;
smth=@(n) conv2(n,[0.05,0.1,0.05;0.1,0.4,0.1;0.05,0.1,0.05;],'same');
% smth=@(n) conv2(n,[0.08,0.15,0.08;0.15,0.23,0.15;0.08,0.15,0.08;],'same');
percents=[0.3,0.6,0.9];
myfigure;
f=gal.CentralSampleIter>0&gal.FlagSFR>0&~isnan(gal.Umodel_abs)&~isnan(gal.Rmodel_abs);
x=gal.color(f);y=log10(gal.SSFR(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[0,5],[-3,1]);
levels=percentile_to_density(n,percents);
[~,h0]=contour(xx,yy,smth(n),levels,'r');hold on;%colormap(cmap);
f=gal.CentralSampleIter>0&~isnan(gal.Umodel_abs)&~isnan(gal.Rmodel_abs);
gal.SSFR2=gal.SSFR;
% gal.SSFR2(gal.FlagSFR<=0)=10.^(-2.6+rand(sum(gal.FlagSFR<=0),1));
gal.SSFR2(gal.SSFR<=0)=nan;%10.^(-2.6+rand(sum(gal.SSFR<=0),1));
x=gal.color(f);y=log10(gal.SSFR2(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[0,5],[-3,1]);
[~,h0]=contour(xx,yy,smth(n),levels,'b--');hold on;%colormap(cmap);
f=galmock.CentralSampleIter>0&galmock.zhub<0.31;
x=galmock.color(f);y=log10(galmock.SSFR(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[0,5],[-3,1]);
[~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'g');hold on;%colormap(cmap);
% f=rawmock.z<0.31;
% x=log10(rawmock.sm(f));y=log10(rawmock.SFR(f)./rawmock.sm(f));
% [xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[-5,4]);
% [~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'b');hold on;%colormap(cmap);
plot([0,5],[-1.5,-1.5],'k--');
plot([2,2],[-4,2],'k--');
plot([1.6,1.6],[-4,2],'k:');
legend('Data','Data: Make-up','Mock');
xlabel('U-R');
ylabel('$\log$(SSFR[Gyr$^{\frac{\ }{\ }1}$])');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/SSFR-Color.eps');
%% SSFR-mag
nbin=40;
smth=@(n) conv2(n,[0.05,0.1,0.05;0.1,0.4,0.1;0.05,0.1,0.05;],'same');
% smth=@(n) conv2(n,[0.08,0.15,0.08;0.15,0.23,0.15;0.08,0.15,0.08;],'same');
percents=[0.3,0.6,0.9];
myfigure;
f=gal.CentralSampleIter>0&gal.FlagSFR>0&~isnan(gal.Umodel_abs)&~isnan(gal.Rmodel_abs);
x=gal.Rmodel_abs(f);y=log10(gal.SSFR(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[-24,-15],[-3,1]);
levels=percentile_to_density(n,percents);
[~,h0]=contour(xx,yy,smth(n),levels,'r');hold on;%colormap(cmap);
f=gal.CentralSampleIter>0&~isnan(gal.Umodel_abs)&~isnan(gal.Rmodel_abs);
gal.SSFR2=gal.SSFR;
gal.SSFR2(gal.FlagSFR<=0)=10.^(-2.6+rand(sum(gal.FlagSFR<=0),1));
x=gal.Rmodel_abs(f);y=log10(gal.SSFR2(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[-24,-15],[-3,1]);
[~,h0]=contour(xx,yy,smth(n),levels,'b--');hold on;%colormap(cmap);
f=galmock.CentralSampleIter>0&galmock.zhub<0.31;
x=galmock.Rmodel_abs(f);y=log10(galmock.SSFR(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[-24,-15],[-3,1]);
[~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'g');hold on;%colormap(cmap);
% f=rawmock.z<0.31;
% x=log10(rawmock.sm(f));y=log10(rawmock.SFR(f)./rawmock.sm(f));
% [xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[-5,4]);
% [~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'b');hold on;%colormap(cmap);
plot([0,5],[-1.5,-1.5],'k--');
plot([2,2],[-4,2],'k--');
plot([1.6,1.6],[-4,2],'k:');
legend('Data','Data: Make-up','Mock');
xlabel('Abs R mag');
ylabel('$\log$(SSFR[Gyr$^{\frac{\ }{\ }1}$])');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/SSFR-Rmag.eps');
%%
f1=gal.CentralSampleIter>0&gal.FlagSFR>0&gal.SSFR>10^-1.5;
f2=gal.CentralSampleIter>0&~(gal.FlagSFR>0&gal.SSFR>10^-1.5);
figure;
y=(gal.Zspec);
f3=(gal.SMsps>1e10&gal.SMsps<2e10);
linhist(y(f1&f3),20,'stairsnorm','b');
hold on;
linhist(y(f2&f3),20,'stairsnorm','r');
%%
nbin=40;
smth=@(n) conv2(n,[0.05,0.1,0.05;0.1,0.4,0.1;0.05,0.1,0.05;],'same');
% smth=@(n) conv2(n,[0.08,0.15,0.08;0.15,0.23,0.15;0.08,0.15,0.08;],'same');
percents=[0.3,0.6,0.9];
myfigure;
f=gal.CentralSampleIter>0&gal.FlagSFR>0&~isnan(gal.Umodel_abs)&~isnan(gal.Rmodel_abs);
x=gal.Zspec(f);y=log10(gal.SSFR(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[0,0.5],[-3,1]);
levels=percentile_to_density(n,percents);
[~,h0]=contour(xx,yy,smth(n),levels,'r');hold on;%colormap(cmap);
f=gal.CentralSampleIter>0&~isnan(gal.Umodel_abs)&~isnan(gal.Rmodel_abs);
gal.SSFR2=gal.SSFR;
gal.SSFR2(gal.FlagSFR<=0)=10.^(-2.6+rand(sum(gal.FlagSFR<=0),1));
x=gal.Zspec(f);y=log10(gal.SSFR2(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[0,0.5],[-3,1]);
[~,h0]=contour(xx,yy,smth(n),levels,'b--');hold on;%colormap(cmap);
f=galmock.CentralSampleIter>0&galmock.zhub<0.31;
x=galmock.Z(f);y=log10(galmock.SSFR(f));
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[0,0.5],[-3,1]);
[~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'g');hold on;%colormap(cmap);
% f=rawmock.z<0.31;
% x=log10(rawmock.sm(f));y=log10(rawmock.SFR(f)./rawmock.sm(f));
% [xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[-5,4]);
% [~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'b');hold on;%colormap(cmap);
plot([0,5],[-1.5,-1.5],'k--');
plot([2,2],[-4,2],'k--');
plot([1.6,1.6],[-4,2],'k:');
legend('Data','Data: Make-up','Mock');
xlabel('z');
ylabel('$\log$(SSFR[Gyr$^{\frac{\ }{\ }1}$])');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/SSFR-Color.eps');
%% Mstar-z
nbin=40;
smth=@(n) conv2(n,[0.05,0.1,0.05;0.1,0.4,0.1;0.05,0.1,0.05;],'same');
% smth=@(n) conv2(n,[0.08,0.15,0.08;0.15,0.23,0.15;0.08,0.15,0.08;],'same');
percents=[0.3,0.6,0.9];
myfigure;
f=gal.CentralSampleIter>0&gal.FlagSFR>0&~isnan(gal.Umodel_abs)&~isnan(gal.Rmodel_abs);
x=log10(gal.SMsps(f));y=gal.Zspec(f);
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[0,0.5]);
levels=percentile_to_density(n,percents);
[~,h0]=contour(xx,yy,smth(n),levels,'r');hold on;%colormap(cmap);
f=galmock.CentralSampleIter>0&galmock.zhub<0.31;
x=log10(galmock.Mstar(f));y=galmock.Z(f);
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[0,0.5]);
[~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'g');hold on;%colormap(cmap);
% f=rawmock.z<0.31;
% x=log10(rawmock.sm(f));y=log10(rawmock.SFR(f)./rawmock.sm(f));
% [xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[-5,4]);
% [~,h1]=contour(xx,yy,smth(n),percentile_to_density(n,percents),'b');hold on;%colormap(cmap);
legend('Data','Mock');
xlabel('$\log(Mstar)$');
ylabel('z');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/Mstar-Z.eps');
%% Mstar-z (SFgalaxy div)
nbin=50;
smth=@(n) conv2(n,[0.05,0.1,0.05;0.1,0.4,0.1;0.05,0.1,0.05;],'same');
percents=[0.3,0.6,0.9];
myfigure;
f=gal.CentralSampleIter>0&gal.IsAGN==1;
x=log10(gal.SMsps(f));y=gal.Zspec(f);
% [xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[0,0.5]);
% % levels=percentile_to_density(n,percents);
% [~,h0]=contour(xx,yy,smth(n),levels,'b');hold on;%colormap(cmap);
plot(x,y,'b.','markersize',1,'color',[0.6,0.6,0.6]);hold on;
f=gal.CentralSampleIter>0&gal.IsAGN==0;
x=log10(gal.SMsps(f));y=gal.Zspec(f);
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[0,0.5]);
levels=percentile_to_density(n,percents);
[~,h0]=contour(xx,yy,smth(n),levels,'b','linewidth',2);hold on;%colormap(cmap);
f=gal.CentralSampleIter>0;%all
x=log10(gal.SMsps(f));y=gal.Zspec(f);
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[0,0.5]);
% levels=percentile_to_density(n,percents);
[~,h0]=contour(xx,yy,smth(n),levels,'r','linewidth',2);hold on;%colormap(cmap);

f=galmock.CentralSampleIter>0;
x=log10(galmock.Mstar(f));y=galmock.Z(f);
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[0,0.5]);
levels=percentile_to_density(n,percents);
[~,h0]=contour(xx,yy,smth(n),levels,'m--','linewidth',2);hold on;%colormap(cmap);

% f=gal.CentralSampleIter>0;
% [xm,ym]=skeleton(log10(gal.SMsps(f)),gal.Zspec(f),30);
% plot(xm,ym,'--');
% plot(xm,0.0228*xm.^2-0.341*xm+1.307,'-');
l=legend('AGN','SF central','All central','Mock central');
set(l,'location','northwest');
xlabel('$\log(Mstar)$');
ylabel('z');
xlim([8,12]);
% print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/Mstar-Z-divSF.eps');
%% Mstar-z (SSFR div)
nbin=40;
smth=@(n) conv2(n,[0.05,0.1,0.05;0.1,0.4,0.1;0.05,0.1,0.05;],'same');
% smth=@(n) conv2(n,[0.08,0.15,0.08;0.15,0.23,0.15;0.08,0.15,0.08;],'same');
percents=[0.3,0.6,0.9];
myfigure;
f=gal.CentralSampleIter>0&gal.IsAGN==0;
%&gal.FlagSFR>1&gal.SSFR>10^-1.5;
x=log10(gal.SMsps(f));y=gal.Zspec(f);
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[0,0.5]);
levels=percentile_to_density(n,percents);
[~,h0]=contour(xx,yy,smth(n),levels,'r');hold on;%colormap(cmap);
f=gal.CentralSampleIter>0;%&~(gal.FlagSFR>1&gal.SSFR>10^-1.5);
x=log10(gal.SMsps(f));y=gal.Zspec(f);
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[0,0.5]);
% levels=percentile_to_density(n,percents);
[~,h0]=contour(xx,yy,smth(n),levels,'g');hold on;%colormap(cmap);

f=galmock.CentralSampleIter>0;%&galmock.SSFR>10^-1.5;
x=log10(galmock.Mstar(f));y=galmock.Z(f);
[xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[0,0.5]);
levels=percentile_to_density(n,percents);
[~,h0]=contour(xx,yy,smth(n),levels,'m--');hold on;%colormap(cmap);
% f=galmock.CentralSampleIter>0&galmock.SSFR<=10^-1.5;
% x=log10(galmock.Mstar(f));y=galmock.Z(f);
% [xx,yy,n]=densitygrid(x,y,[nbin,nbin],[8,12],[0,0.5]);
% % levels=percentile_to_density(n,percents);
% [~,h0]=contour(xx,yy,smth(n),levels,'c--');hold on;%colormap(cmap);

f=gal.CentralSampleIter>0;
[xm,ym]=skeleton(log10(gal.SMsps(f)),gal.Zspec(f),30);
plot(xm,ym,'--');
plot(xm,0.0228*xm.^2-0.341*xm+1.307,'-');
l=legend('Active','Inactive','Mock Active','Mock Inactive');
set(l,'location','northwest');
xlabel('$\log(Mstar)$');
ylabel('z');
xlim([8,12]);
% print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/Mstar-Z-div.eps');
%% 
myfigure;
fz=gal.Zspec>0.0228*log10(gal.SMsps).^2-0.341*log10(gal.SMsps)+1.307;
[xm1,ym1]=linhist(log10(gal.SMsps(fz&gal.CentralSampleIter>0&(gal.FlagSFR>0&gal.SSFR>10^-1.5))),[xbrd(:,1);xbrd(end)])
[xm,ym]=linhist(log10(gal.SMsps(fz&gal.CentralSampleIter>0)),[xbrd(:,1);xbrd(end)]);
plot(xm,ym1./ym,'g-o');hold on;
[xm1,ym1]=linhist(log10(gal.SMsps(gal.CentralSampleIter>0&(gal.FlagSFR>0&gal.SSFR>10^-1.5))),[xbrd(:,1);xbrd(end)])
[xm,ym]=linhist(log10(gal.SMsps(gal.CentralSampleIter>0)),[xbrd(:,1);xbrd(end)]);
plot(xm,ym1./ym,'r-o');hold on;
% yscale('log');
xlabel('$M_\star[\rm{M}_\odot/h^2]$');
ylabel('Active Fraction');
legend('Highz sample','whole sample');
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/ActiveFraction-Zsample.eps');
%% compare mock group SSFR dependence
myfigure;
f=grpmock.SSFRIter>10^-1.5; %about 90% are real centrals
x=logspace(8,12,20);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.MstarIter(f),log10(grpmock.Luminosity(f)),x,0.68);
h2=plot(xmed,ymed,'--','linewidth',2,'color','g');hold on;
f=grpmock.Mult>0;
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.MstarIter(f),log10(grpmock.Luminosity(f)),x,0.68);
h2=plot(xmed,ymed,'--','linewidth',2,'color','r');hold on;

f=grp.IterCenSSFR>10^-1.5; %about 90% are real centrals
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grp.IterCenSM(f),log10(grp.Luminosity(f)),x,0.68);
h2=plot(xmed,ymed,'-','linewidth',2,'color','g');hold on;
f=grp.Mult>0;
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grp.IterCenSM(f),log10(grp.Luminosity(f)),x,0.68);
h2=plot(xmed,ymed,'-','linewidth',2,'color','r');hold on;
xscale('log');
legend('Mock Active','Mock All','data active','data all');
xlim([5e9,1e12]);
%% rmag dependence
f=bower.sm>1e10&bower.sm<2e10;
mbrd=(minmax(bower.Rmag(f)'));
nbin=linspace(mbrd(1),mbrd(2),40);
colors=jet(numel(nbin));
myfigure;
for i=1:numel(nbin)-1
    ff=f&bower.Rmag>=nbin(i)&bower.Rmag<nbin(i+1);
    loglog(bower.sm(ff),bower.mh(ff),'o','color',colors(i,:));
hold on;
end
xscale('linear');
yscale('log');
xlabel('$M_\star [\rm{M}_\odot/h]$');
ylabel('$M_{halo} [\rm{M}_\odot/h]$');
% xlim([0.85e11,1.05e11]);
colormap(jet(numel(nbin)));%caxis([1, numel(nbin)]);
hcb = colorbar('YTick',(([-22.1,-21,-20,-19])-mbrd(1))/(mbrd(2)-mbrd(1))*(numel(nbin)-1)+1,'YTickLabel',...
 {'Mag','-21','-20','-19'});
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/Mhalo-Mstar-Rmag.eps');
%%
myfigure;
f=bower.sm>1e10&bower.sm<1.3e10;
loglog(bower.Rmag(f),bower.mh(f),'.');
xscale('linear');
yscale('log');
xlabel('$M_{r}$');
ylabel('$M_{halo} [\rm{M}_\odot/h]$');
% xlim([0.85e11,1.05e11]);
% colormap(jet(numel(nbin)));%caxis([1, numel(nbin)]);
% hcb = colorbar('YTick',(log10([1,10,100,1000])-mbrd(1))/(mbrd(2)-mbrd(1))*(numel(nbin)-1)+1,'YTickLabel',...
% {'N','10','100','1000'});
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/Mhalo-Rmag-FixMstar.eps');
%%
%% compare with grp measurement
dataZ02=[%7.28742e+09   1.45184e+09  11.63      0.5295                  
1.43904e+10   2.89589e+09  11.72      0.5876                  
2.52734e+10   3.55315e+09  12.16      0.3631                  
4.10653e+10   5.77329e+09  12.47      0.2683                  
6.55599e+10   9.16918e+09  12.79      0.2516                  
1.05014e+11   1.43846e+10  13.65      0.1315                  
1.67394e+11   2.12134e+10  13.6       0.2977                  
2.74124e+11   2.95823e+10  13.33      1.34                    
5.09854e+11   7.66566e+10  8.618      nan
]; %z<0.2
xbrd=[1e10,2e10; 0.2e11,0.325e11; 0.3251e11,0.5285e11; 0.5285e11,0.8592e11; 0.8592e11,1.3967e11; 1.3967e11,2.2705e11; 2.2705e11,3.7e11;3.7e11, 6e11];
% xbrd=log10(xbrd);
myfigure;
sigma_common=0.2;
for i=1:numel(ms)
    zwang(i)=medianMh(log10(ms(i)),sigma_common,Pwang4);
    zling(i)=medianMh(log10(ms(i)),sigma_common,Pling);
    zmosterZ0(i)=medianMh(log10(ms(i)),sigma_common,PmosterZ0);
    zguo(i)=medianMh(log10(ms(i)),sigma_common,Pguo); 
end
zbound=minmax([zwang',zling',zmosterZ0',zguo']);
h3=area(ms,[zbound(:,1),zbound(:,2)-zbound(:,1)]);hold on;
set(h3(2),'facecolor','y','edgecolor','w');
set(h3(1),'facecolor','w','edgecolor','w');
set(gca,'layer','top');

% f=galmock.CentralSampleIter>0&galmock.Z<0.2;
% [xmed,ymed,yl,xm,ym,yerr]=skeleton(galmock.Mstar(f),log10(galmock.Mhalo(f)),x,0.68);
% h11=plot(xmed,ymed,'k--');
h4=plot(ms,ywang4,'c--');
h5=plot(ms,yling,'m--');
h6=plot(ms,ymosterZ0,'g--');
mh=logspace(10,16,10);
h7=plot(halo2starP(mh,Pguo),log10(mh),'b--'); 

h1=ploterr(dataZ02(:,1),dataZ02(:,3)-0.,{xbrd(:,1),xbrd(:,2)},dataZ02(:,4),'ks','logx');hold on;
xlim([9e9,1e12]);
ylim([10.8,15.5]);
xscale('log');
xlabel('$M_\star$[M$_\odot/h^2$]');
ylabel('$\log(M_h$[M$_\odot/h$])');
l=legend([h1(1),h4,h5,h6,h7,h3(2)],'GAMA (this work)','WangLan13','WangLingyu13','Moster13','Guo10','All Models $\sigma=0.2$');
% h4=plot(ms,ywang4,'c--');
% h5=plot(ms,yling,'m--');
% h6=plot(ms,ymosterZ0,'g--');
% mh=logspace(10,16,10);
% h7=plot(halo2starP(mh,Pguo),log10(mh),'b--'); 

set(l,'location','northwest','interpreter','latex');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/extra/M-Mstar-poster.eps');