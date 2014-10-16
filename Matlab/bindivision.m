% clear;
% clear global;
load_G3Cv4;

msun=4.67; %solar ab abs mag
% L/Lsun=10^(-0.4*(m-msun));
% grp.LBCG=10.^(-0.4*(grp.BCGmag-msun));
% grp.L=grp.L*0.7^2;
% grp.LBCG=grp.LBCG*0.7^2;

bwl=[-0.01,-0.04,-0.09];  %WL photo-z calibration bias. signal_true=signal_obs/(1+bwl).
%%
f=grp.Mult>2;
hybmass=grp.HybMass(f);
dynmass=grp.MassProxy(f);
lummass=grp.LumMass(f);
veldisp=grp.VelDisp(f);
z=grp.Zfof(f);
L=grp.TotFluxProxy(f);
% S=grp.TotSMInt(f);
% G=grp.LBCG(f);
m=lummass;
[m,ix]=sort(m);
z=z(ix);
hybmass=hybmass(ix);
dynmass=dynmass(ix);
lummass=lummass(ix);
L=L(ix);
% S=S(ix);
% G=G(ix);

n=numel(m);
n1=floor(n/3);n2=floor(2*n/3);

md=zeros(2,2);
md(1,:)=(m(floor(numel(m)/3)+(0:1)));
md(2,:)=(m(floor(numel(m)*2/3)+(0:1)));
md

% figure;
% nbin=31;
% [xm1,ym1]=linhist(z(1:n1),nbin,'stairs','r');
% hold on;
% [xm2,ym2]=linhist(z(n1+1:n2),nbin,'stairs','g');
% [xm3,ym3]=linhist(z(n2+1:n),nbin,'stairs','b');
% dat=[xm1,ym1,xm2, ym2,xm3, ym3];
% dlmwrite('LumBinZdistr.dat',dat);
% dlmwrite('VelDispBinZdistr.dat',dat);
% dlmwrite('DynMassZdistr.dat',dat);
% dlmwrite('LumMassBinZdistr.dat',dat);

% dlmwrite('L3xyz.dat',[grp.BCGRA(n2+1:n),grp.BCGDEC(n2+1:n),grp.Zfof(n2+1:n)]);

% fun=@mean; 
% efun=@(x) fun(x)+[-std(x)/sqrt(numel(x)); std(x)/sqrt(numel(x))]; %error on the mean
fun=@median;
efun=@(x) (prctile(x,[(100-68.3)/2;100-(100-68.3)/2])-fun(x))/sqrt(numel(x))+fun(x);
Zbin=[fun(z(1:n1)),fun(z(n1+1:n2)),fun(z(n2+1:n))];
Lbin=[fun(L(1:n1)),fun(L(n1+1:n2)),fun(L(n2+1:n))];
% Sbin=[fun(S(1:n1)),fun(S(n1+1:n2)),fun(S(n2+1:n))];
% Gbin=[fun(G(1:n1)),fun(G(n1+1:n2)),fun(G(n2+1:n))];
%m=hybmass;
%f=grp.Mult>2;m1=m(1:n1);m2=m(n1+1:n2);m3=m(n2+1:n);
%m1=m1(f(1:n1));m2=m2(f(n1+1:n2));m3=m3(f(n2+1:n));
%MHbin=[fun(m1),fun(m2),fun(m3)];
% fun=@median;efun=@(x) prctile(x,[(100-68.3)/2;100-(100-68.3)/2]);
m=hybmass;
MHbin=[fun(m(1:n1)),fun(m(n1+1:n2)),fun(m(n2+1:n))];
eMHbin=[efun(m(1:n1)),efun(m(n1+1:n2)),efun(m(n2+1:n))];
m=dynmass;
MDbin=[fun(m(1:n1)),fun(m(n1+1:n2)),fun(m(n2+1:n))];
eMDbin=[efun(m(1:n1)),efun(m(n1+1:n2)),efun(m(n2+1:n))];
m=lummass;
MLbin=[fun(m(1:n1)),fun(m(n1+1:n2)),fun(m(n2+1:n))];
eMLbin=[efun(m(1:n1)),efun(m(n1+1:n2)),efun(m(n2+1:n))];
eLbin=[efun(L(1:n1)),efun(L(n1+1:n2)),efun(L(n2+1:n))];
ML_L=MLbin./Lbin;
eML_L=[ML_L;ML_L]+sqrt((eMLbin./[MLbin;MLbin]-1).^2+(eLbin./[Lbin;Lbin]-1).^2).*[-ML_L;ML_L];
MD_L=MDbin./Lbin;
eMD_L=[MD_L;MD_L]+sqrt((eMDbin./[MDbin;MDbin]-1).^2+(eLbin./[Lbin;Lbin]-1).^2).*[-MD_L;MD_L];
figure;
loglog(lummass,dynmass,'.');
hold on;
ploterr(MLbin,MDbin,{eMLbin(1,:),eMLbin(2,:)},{eMDbin(1,:),eMDbin(2,:)},'ro');
% ploterr(MLbin,MDbin,eMLbin,eMDbin,'ro');
plot([1e12,1e15],[1e12,1e15],'r')


%%
sum(m<md(1))/numel(m)
sum(m>md(1)&m<md(2))/numel(m)
sum(m>md(2))/numel(m)
%%
OmegaM=0.25;H0=100;G=43.0071;
rhob=OmegaM*3*H0^2/8/pi/G;
binname='LumMass';

myfigure;
[h,r1,s1,es1,cov1]=plotgama_prof('L1V4','r.');
hold on;
[h,r2,s2,es2,cov2]=plotgama_prof('L2V4','go');
[h,r3,s3,es3,cov3]=plotgama_prof('L3V4','bs');
%%
s1=s1/(1+bwl(1)); es1=es1/(1+bwl(1));
s2=s2/(1+bwl(2)); es2=es2/(1+bwl(2));
s3=s3/(1+bwl(3)); es3=es3/(1+bwl(3));
%%
fitfun=@createFit_auto_concentration;
% fitfun=@createFit;
[cf1,gof1,out1]=fitfun(r1,s1,es1.^-2,Zbin(1));
ci1=confint(cf1,0.683);
par1=coeffvalues(cf1);
[cf2,gof2,out2]=fitfun(r2,s2,es2.^-2,Zbin(2));
ci2=confint(cf2,0.683);
par2=coeffvalues(cf2);
[cf3,gof3,out3]=fitfun(r3,s3,es3.^-2,Zbin(3));
ci3=confint(cf3,0.683);
par3=coeffvalues(cf3);

mfit=[par1(1);par2(1);par3(1)]*1e10;
bfit=[par1(2);par2(2);par3(2)];
% cfit=[par1(3);par2(3);par3(3)];
errm=[ci1(:,1),ci2(:,1),ci3(:,1)]*1e10;
errb=[ci1(:,2),ci2(:,2),ci3(:,2)];
% errc=[ci1(:,3),ci2(:,3),ci3(:,3)];
%%
fitfun=@createFit;
[cf1,gof1,out1]=fitfun(r1,s1,es1.^-2,Zbin(1),MLbin(1)/1e10);
ci1=confint(cf1,0.683);
par1=coeffvalues(cf1);
[cf2,gof2,out2]=fitfun(r2,s2,es2.^-2,Zbin(2),MLbin(2)/1e10);
ci2=confint(cf2,0.683);
par2=coeffvalues(cf2);
[cf3,gof3,out3]=fitfun(r3,s3,es3.^-2,Zbin(3),MLbin(3)/1e10);
ci3=confint(cf3,0.683);
par3=coeffvalues(cf3);

bfit=[par1(1);par2(1);par3(1)];
cfit=[par1(2);par2(2);par3(2)];
errb=[ci1(:,1),ci2(:,1),ci3(:,1)];
errc=[ci1(:,2),ci2(:,2),ci3(:,2)];
%%
J=out1.Jacobian;
ok= isfinite(r1) & isfinite(s1) & isfinite(es1);  % to account for number of observations in fitting
err1=J\cov1(ok,ok)/J';
err=inv(J'*J);
J=out2.Jacobian;
ok= isfinite(r2) & isfinite(s2) & isfinite(es2);
err2=J\cov2(ok,ok)/J';
J=out3.Jacobian;
ok= isfinite(r3) & isfinite(s3) & isfinite(es3);
err3=J\cov3(ok,ok)/J';
em=sqrt([err1(1,1),err2(1,1),err3(1,1)])';
eb=sqrt([err1(2,2),err2(2,2),err3(2,2)])';
ec=sqrt([err1(3,3),err2(3,3),err3(3,3)])';
cmb=[err1(1,2),err2(1,2),err3(1,2)]';  %correlation between m and b
cmb=cmb./em./eb;
% plot(r1,log10(nfw_DeltSig(r1./(1+Zbin(1)),Mbin(1)/1e10,Zbin(1),2)/(1+Zbin(1))^2),'r-');
% plot(r2,log10(nfw_DeltSig(r2./(1+Zbin(2)),Mbin(2)/1e10,Zbin(2),2)/(1+Zbin(2))^2),'g-');
% plot(r3,log10(nfw_DeltSig(r3./(1+Zbin(3)),Mbin(3)/1e10,Zbin(3),2)/(1+Zbin(3))^2),'b-');
% plot(r1,log10(2*rhob*ppval(DSig_interp_spline,r1).*growth_factor(0.3,Zbin(1)).^2),'r--');
% plot(r2,log10(2.5*rhob*ppval(DSig_interp_spline,r2).*growth_factor(0.3,Zbin(2)).^2),'g--');
% plot(r3,log10(4*rhob*ppval(DSig_interp_spline,r3).*growth_factor(0.3,Zbin(3)).^2),'b--');

plot(r1,signlog(lensing_rfunc(r1,par1(1),par1(2),Zbin(1),par1(3))),'r-'); hold on;
plot(r2,signlog(lensing_rfunc(r2,par2(1),par2(2),Zbin(2),par2(3))),'g-');
plot(r3,signlog(lensing_rfunc(r3,par3(1),par3(2),Zbin(3),par3(3))),'b-');
%%
myfigure;
set(gcf,'defaultlinemarkersize',8);
l1=ploterr(r1,s1,[],es1,'rd','logxy');hold on;set(l1(1),'markerfacecolor','r');
ltmp=ploterr(r1(s1<0),-s1(s1<0),[],[],'rd','logxy');%set(ltmp(1),'markerfacecolor','r');
l2=ploterr(r2,s2,[],es2,'go','logxy');set(l2(1),'markerfacecolor','g');
ltmp=ploterr(r2(s2<0),-s2(s2<0),[],[],'go','logxy');%set(ltmp(1),'markerfacecolor','b');
l3=ploterr(r3,s3,[],es3,'bs','logxy');set(l3(1),'markerfacecolor','b');
ltmp=ploterr(r3(s3<0),-s3(s3<0),[],[],'bs','logxy','abshhy',0.1);%set(ltmp(1),'markerfacecolor','m');
% plot(r1,(nfw_DeltSig(r1./(1+Zbin(1)),Mbin(1)/1e10,Zbin(1),2)/(1+Zbin(1))^2),'r-');
% plot(r2,(nfw_DeltSig(r2./(1+Zbin(2)),Mbin(2)/1e10,Zbin(2),2)/(1+Zbin(2))^2),'b-');
% plot(r3,(nfw_DeltSig(r3./(1+Zbin(3)),Mbin(3)/1e10,Zbin(3),2)/(1+Zbin(3))^2),'m-');
% plot(r1,(lensing_rfunc(r1,par1(1),par1(2),Zbin(1),par1(3))),'r-'); hold on;
% plot(r2,(lensing_rfunc(r2,par2(1),par2(2),Zbin(2),par2(3))),'g-');
% plot(r3,(lensing_rfunc(r3,par3(1),par3(2),Zbin(3),par3(3))),'b-');
plot(r1,(lensing_rfunc(r1,mfit(1)/1e10,bfit(1),Zbin(1),cfit(1))),'r-'); hold on;
plot(r2,(lensing_rfunc(r2,mfit(2)/1e10,bfit(2),Zbin(2),cfit(2))),'g-');
plot(r3,(lensing_rfunc(r3,mfit(3)/1e10,bfit(3),Zbin(3),cfit(3))),'b-');
% set(gca,'xscale','log','yscale','log');
makeup_errbar(r1,s1,es1,'r-');
makeup_errbar(r2,s2,es2,'g-');
makeup_errbar(r3,s3,es3,'b-');
l=legend([l1(1),l2(1),l3(1)],...
    [printexp10(mfit(1),'%1.0f'),', ',num2str(Zbin(1),'%.2f')],...
    [printexp10(mfit(2),'%1.0f'),', ',num2str(Zbin(2),'%.2f')],...
    [printexp10(mfit(3),'%1.0f'),', ',num2str(Zbin(3),'%.2f')]);
% set(l,'box','off');
xlabel('r/(Mpc/h) comoving');
ylabel('$\Delta\Sigma/(10^{10}M_\odot/h/(Mpc/h)^2)$');
% print('-depsc',['wlSig_fit_',binname,'.eps']);
%%
h=1;
myfigure;set(gcf,'defaultlinemarkersize',8);
% HH=ploterr(mfit,MHbin*h,{errm(1,:)',errm(2,:)'},{eMHbin(1,:)*h,eMHbin(2,:)*h},'o','logxy');hold on;
HL=ploterr(mfit,MLbin*h,em*1e10,{eMLbin(1,:)*h,eMLbin(2,:)*h},'r>','logxy');hold on;
HD=ploterr(mfit,MDbin*h,em*1e10,{eMDbin(1,:)*h,eMDbin(2,:)*h},'gs','logxy');
% HS=ploterr(mfit,Sbin,{errm(1,:)',errm(2,:)'},[],'k*','logxy');
l=legend([HL(1),HD(1)],'LumMass','DynMass','SteMass');
set(l,'location','northwest');
plot([1e11,1e15], [1e11,1e15],'--');
xlabel('$M_{WL}/(M_\odot/h)$');ylabel('$M_{est}/(M_\odot/h)$, Median');
xlim([1e11,1e15]);ylim([1e11,1e15]);
set(gca,'xtick',logspace(11,15,5));
% print('-depsc',['MM_',binname,'median.eps']);
%% compare mass concentration
% allowing concentration free doesn't solve the problem; and the M-c
% relation is consistent with theoretical predictions
mfit=MLbin';em=0;
l1=ploterr(mfit,cfit,em*1e10,ec,'o','logx','logy');
hold on;
Mc=2e12
cguess=10.14.*(1+Zbin).^-1.01.*(mfit'/Mc).^-0.081
l2=plot(mfit,cguess,'-')
x=logspace(12,14,5);
l3=plot(x,10.14.*(x/Mc).^-0.081,'r');
legend([l1(1),l2,l3],'fit','Duffy with z','Duffy');
%%
myfigure;set(gcf,'defaultlinemarkersize',8);
sig8=0.6;
x=logspace(11,14,5);
y1=linear_bias(x/1e10,0.2,0.6);
y2=linear_bias(x/1e10,0.2,0.8);
y3=linear_bias(x/1e10,0.2,1);
% l=ploterr(mfit,bfit,{errm(1,:)',errm(2,:)'},{errb(1,:)',errb(2,:)'},'o','logx');
% hold on;
l=ploterr(mfit,bfit,em*1e10,eb,'ro','logx');hold on;
% l=ploterr(MLbin,bfit,{eMLbin(1,:)*h,eMLbin(2,:)*h},eb,'gs','logx');hold on;
lt=[];
h=plot(x,y1,'k:','displayname','\sigma_8=0.6');    lt=[lt;h(1)];
h=plot(x,y2,'k-.','displayname','\sigma_8=0.8');    lt=[lt;h(1)];
h=plot(x,y3,'k-','displayname','\sigma_8=1');    lt=[lt;h(1)];
legend(lt,'location','southwest');
xlabel('$M_{WL}$');ylabel('bias');%legend(l(1),'HybMass bin','location','northwest');
print('-depsc',['bias_',binname,'.eps']);
%%
myfigure;
ml=[9.2,2.265;9.5,2.2;9.8,2.08;10.1,2.035;10.4,1.95;10.7,2.07;11,2.22;11.3,2.4;11.6,2.48;11.88,2.55];
ml(:,2)=ml(:,2);
ml=10.^ml;
% h2=ploterr(Lbin,MD_L,eLbin,{eMDbin(1,:)./Lbin,eMDbin(2,:)./Lbin},'gs','logxy');hold on;
h2=ploterr(Lbin,MD_L,{eLbin(1,:),eLbin(2,:)},{eMD_L(1,:),eMD_L(2,:)},'gs','logxy');hold on;
h1=ploterr(Lbin,mfit./Lbin',{eLbin(1,:),eLbin(2,:)},em*1e10./Lbin','ro','logxy');
hold on;

h3=plot(ml(:,1),ml(:,2),'s-');
xlabel('$L_r/(h^{\frac{\ }{\ }2}L_{r,{\odot}})$');ylabel('$M/L_r [h M_\odot/L_{r,\odot}]$');legend([h1(1),h2(1),h3],'WL','DynMass','2PIGG','location','northwest');
print('-depsc',['ML_',binname,'.eps']);
%%
myfigure;
% h1=ploterr(mfit,100*Gbin,{errm(1,:)',errm(2,:)'},[],'ro','logxy');
% hold on;
h2=ploterr(mfit,100*Lbin,{errm(1,:)',errm(2,:)'},[],'gs','logxy');hold on;
xlabel('$M_{WL}$');ylabel('$100*L$');
plot([1e12,1e14], [1e12,1e14],'--');
% legend([h1(1),h2(1)],'BCG','Total','location','northwest');
% print('-depsc',['ML2_',binname,'.eps']);
%%
% allfits=[];

savefit.name=binname;
savefit.cf={cf1;cf2;cf3};
savefit.gof={gof1;gof2;gof3};
savefit.mfit=mfit;
savefit.bfit=bfit;
savefit.errm=errm';
savefit.errb=errb';
savefit.mhyb=MHbin';
allfits=[allfits;savefit];

%  save allfits.mat allfits
%%
%load allfits.mat
linespec={'ro','gs','b>'};
myfigure;
l=[];
for i=1:numel(allfits)
    mh=allfits(i).mhyb;
    m=allfits(i).mfit;
    em=allfits(i).errm;
    b=allfits(i).bfit;
    eb=allfits(i).errb;
    h=ploterr(m,mh,{em(:,1),em(:,2)},[],linespec{i},'logxy');
    hold on;
    set(h(1),'displayname',allfits(i).name);
    l=[l;h(1)];
end
legend(l,'location','northwest');
plot([1e11,1e14], [1e11,1e14],'--');
xlabel('$M_{WL}$');ylabel('$M_{hyb}$');
% print('-depsc',['MM_allbins.eps']);
%%
myfigure;
l=[];
sig8=0.6;
x=logspace(11,14,5);
y1=linear_bias(x/1e10,0.2,0.6);
y2=linear_bias(x/1e10,0.2,0.8);
y3=linear_bias(x/1e10,0.2,1);
for i=1:numel(allfits)
    mh=allfits(i).mhyb;
    m=allfits(i).mfit;
    em=allfits(i).errm;
    b=allfits(i).bfit;
    eb=allfits(i).errb;
    h=ploterr(m,b,{em(:,1),em(:,2)},{eb(:,1),eb(:,2)},linespec{i},'logx');
    hold on;
    set(h(1),'displayname',allfits(i).name);
    l=[l;h(1)];
end
h=plot(x,y1,'k:','displayname','\sigma_8=0.6');    l=[l;h(1)];
h=plot(x,y2,'k-.','displayname','\sigma_8=0.8');    l=[l;h(1)];
h=plot(x,y3,'k-','displayname','\sigma_8=1');    l=[l;h(1)];
legend(l,'location','southwest');
xlabel('$M_{WL}$');ylabel('$b_{gm}$');
% print('-depsc',['bias_allbins.eps']);

%%
myfigure;
errorbar(r1,s1,es1,'r.');
hold on;
errorbar(r2,s2,es2,'bo');
errorbar(r3(s3>0),s3(s3>0),es3(s3>0),'ms');
errorbar(r3(s3<0),-s3(s3<0),es3(s3<0),'ms','markerfacecolor','m');
plot(r1,(nfw_DeltSig(r1./(1+Zbin(1)),Mbin(1)/1e10,Zbin(1),2)/(1+Zbin(1))^2),'r-');
plot(r2,(nfw_DeltSig(r2./(1+Zbin(2)),Mbin(2)/1e10,Zbin(2),2)/(1+Zbin(2))^2),'b-');
plot(r3,(nfw_DeltSig(r3./(1+Zbin(3)),Mbin(3)/1e10,Zbin(3),2)/(1+Zbin(3))^2),'m-');
set(gca,'xscale','log','yscale','log');
makeup_errbar(r1,s1,es1,'r-');
makeup_errbar(r2,s2,es2,'b-');
makeup_errbar(r3(s3>0),s3(s3>0),es3(s3>0),'m-');
makeup_errbar(r3(s3<0),s3(s3<0),es3(s3<0),'m--');
%%
myfigure;
l1=ploterr(r1,s1,[],es1,'r*','logxy');
hold on;
l2=ploterr(r2,s2,[],es2,'bo','logxy');
l3=ploterr(r3(s3>0),s3(s3>0),[],es3(s3>0),'ms','logxy');
l4=ploterr(r3(s3<0),-s3(s3<0),[],es3(s3<0),'ms','logxy','abshhy',0.1);set(l4(1),'markerfacecolor','m');
plot(r1,(nfw_DeltSig(r1./(1+Zbin(1)),Mbin(1)/1e10,Zbin(1),2)/(1+Zbin(1))^2),'r-');
plot(r2,(nfw_DeltSig(r2./(1+Zbin(2)),Mbin(2)/1e10,Zbin(2),2)/(1+Zbin(2))^2),'b-');
plot(r3,(nfw_DeltSig(r3./(1+Zbin(3)),Mbin(3)/1e10,Zbin(3),2)/(1+Zbin(3))^2),'m-');
% set(gca,'xscale','log','yscale','log');
makeup_errbar(r1,s1,es1,'r-');
makeup_errbar(r2,s2,es2,'b-');
makeup_errbar(r3(s3>0),s3(s3>0),es3(s3>0),'m-');
makeup_errbar(r3(s3<0),-s3(s3<0),es3(s3<0),'m--');
l=legend([l1(1),l2(1),l3(1)],...
    [printexp10(Mbin(1),'%1.0f'),', ',num2str(Zbin(1),'%.2f')],...
    [printexp10(Mbin(2),'%1.0f'),', ',num2str(Zbin(2),'%.2f')],...
    [printexp10(Mbin(3),'%1.0f'),', ',num2str(Zbin(3),'%.2f')]);
% set(l,'box','off');
xlabel('r/(Mpc/h) comoving');
ylabel('$\Delta\Sigma/(10^{10}M_\odot/h/(Mpc/h)^2)$');
% print('-depsc','wlSig_nfwcomp.eps');
%%
grpmock.HybMass=sqrt(grpmock.LumMass.*grpmock.MassProxy);
f=grpmock.Mult>2;%&mockmass.MIter>1e9;
halomass=mockmass.MIter(f);
hybmass=grpmock.HybMass(f);
dynmass=grpmock.MassProxy(f);
lummass=grpmock.LumMass(f);
veldisp=grpmock.VelDisp(f);
z=grpmock.Zfof(f);
L=grpmock.TotFluxProxy(f);
% S=grp.TotSMInt(f);
% G=grp.LBCG(f);
m=L;
[m,ix]=sort(m);
z=z(ix);
halomass=halomass(ix);
hybmass=hybmass(ix);
dynmass=dynmass(ix);
lummass=lummass(ix);
L=L(ix);
% S=S(ix);
% G=G(ix);

n=numel(m);
n1=floor(n/3);n2=floor(2*n/3);

md=zeros(2,2);
md(1,:)=(m(floor(numel(m)/3)+(0:1)));
md(2,:)=(m(floor(numel(m)*2/3)+(0:1)));
md

fun=@median;

Zbin=[fun(z(1:n1)),fun(z(n1+1:n2)),fun(z(n2+1:n))];
Lbin=[fun(L(1:n1)),fun(L(n1+1:n2)),fun(L(n2+1:n))];
% Sbin=[fun(S(1:n1)),fun(S(n1+1:n2)),fun(S(n2+1:n))];
% Gbin=[fun(G(1:n1)),fun(G(n1+1:n2)),fun(G(n2+1:n))];
%m=hybmass;
%f=grp.Mult>2;m1=m(1:n1);m2=m(n1+1:n2);m3=m(n2+1:n);
%m1=m1(f(1:n1));m2=m2(f(n1+1:n2));m3=m3(f(n2+1:n));
%MHbin=[fun(m1),fun(m2),fun(m3)];
m=hybmass;
MHbin=[fun(m(1:n1)),fun(m(n1+1:n2)),fun(m(n2+1:n))];
m=dynmass;
MDbin=[fun(m(1:n1)),fun(m(n1+1:n2)),fun(m(n2+1:n))];
m=lummass;
MLbin=[fun(m(1:n1)),fun(m(n1+1:n2)),fun(m(n2+1:n))];
m=halomass;
Mbin=[fun(m(1:n1)),fun(m(n1+1:n2)),fun(m(n2+1:n))];

% myfigure;
loglog(Mbin,MLbin,'ro');%,'markerfacecolor','r');
hold on;
% loglog(Mbin,MHbin,'x');
loglog(Mbin,MDbin,'bs');%,'markerfacecolor','b');
xlabel('$M_{halo}$');ylabel('$M_{proxy}$');legend('LumMass All','DynMass All','LumMass','DynMass');
plot([1e12,1e15],[1e12,1e15])
title('Mult>1');
% print('-depsc','MassCalibration1.eps');