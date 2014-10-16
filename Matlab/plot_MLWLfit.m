cd /work/Projects/Lensing/outputv4/paper
%% raw fit
p0=3.7;
dynA0=p0*grp.MassProxyRaw;
p1=[8.4,-16.1,0.9];
dynF0=(p1(1)+p1(2)./sqrt(grp.Mult)+p1(3)./sqrt(grp.Zfof)).*grp.MassProxyRaw;
p2=[0.05,-0.27];
dynB0=10^p2(1)*(grp.MassProxy/1e12).^p2(2).*grp.MassProxy;
p3=[-0.29,0.03];
% p3=[-0.14,-0.03];
lumB0=10^p3(1)*(grp.LumMass/1e12).^p3(2).*grp.LumMass;
p4=[14.3,1.11];
lumM0=10^p4(1)*(grp.TotFluxProxy/1e12).^p4(2);
p5=[15.2,1.64];
starM0=10^p5(1)*(grp.IterCenSM/1e12).^p5(2);
p3N5=[-0.14,-0.03];
lumB0N5=10^p3N5(1)*(grp.LumMass/1e12).^p3N5(2).*grp.LumMass;
p3N10=[log10(0.989629) -0.008871];
lumB0N10=10^p3N10(1)*(grp.LumMass/1e12).^p3N10(2).*grp.LumMass;
%% bias correction
p0=3.7-0.4;
dynA=p0*grp.MassProxyRaw;
p1=[8.4+0.2,-16.1+0.5,0.9-0.4];
dynF=(p1(1)+p1(2)./sqrt(grp.Mult)+p1(3)./sqrt(grp.Zfof)).*grp.MassProxyRaw;
p2=[0.05-0.12,-0.27+0.04];
dynB=10^p2(1)*(grp.MassProxy/1e12).^p2(2).*grp.MassProxy;
p3=[-0.29-0.09,0.03+0.01];
lumB=10^p3(1)*(grp.LumMass/1e12).^p3(2).*grp.LumMass;
p4=[14.3-0,1.11+0.01];
lumM=10^p4(1)*(grp.TotFluxProxy/1e12).^p4(2);
p5=[15.2-0.1,1.64+0.01];
starM=10^p5(1)*(grp.IterCenSM/1e12).^p5(2);
%%
f=grp.Mult>2;
m=grp.LumMass(f);
loglog(m,dynA(f),'r.');
hold on;
plot(m,dynF(f),'g.');
plot(m,dynB(f),'b.');
plot(m,lumB(f),'c.');
plot(m,lumM(f),'m.');
plot(m,starM(f),'k.');
%%
figure;
f=grp.Mult>2&grp.Mult<5;
loglog(grp.MassProxy(f),dynF(f),'r.'); hold on;
f=grp.Mult>=5&grp.Mult<10;
loglog(grp.MassProxy(f),dynF(f),'g.')
f=grp.Mult>=10&grp.Mult<20;
loglog(grp.MassProxy(f),dynF(f),'b.')
f=grp.Mult>=20;
loglog(grp.MassProxy(f),dynF(f),'k.')
hold on;
h=legend('[3,5)','[5,10)','[10,20)','[20,1000)')
set(h,'location','southeast');
plot([1e10,1e16],[1e10,1e16])
xlabel('G3CDynMass');
ylabel('Lensing Calibrated Mass');
print('-depsc','DynFreeMass.eps');
%%
figure;
x=logspace(10,16,20);
loghist(grp.MassProxy(f),x,'line','r-');
hold on;
loghist(dynA0(f),x,'line','r-x');
loghist(dynF0(f),x,'line','r-o');
loghist(dynB0(f),x,'line','r--');
loghist(grp.LumMass(f),x,'line','g-');
loghist(lumM0(f),x,'line','g-o');
loghist(lumB0(f),x,'line','g--');
loghist(starM0(f),x,'line','b-');
legend('DynMass','DynA','DynFree','DynBias','LumMass','LumML','LumBias','Star');
xlabel('Mass [Msun/h]');
ylabel('Count');
print('-depsc','MLWLMassDistribution.eps');
%%
figure;
x=logspace(10,16,20);
loghist(grp.MassProxy(f),x,'line','r-');
hold on;
loghist(dynA(f),x,'line','r-x');
loghist(dynF(f),x,'line','r-o');
loghist(dynB(f),x,'line','r--');
loghist(grp.LumMass(f),x,'line','g-');
loghist(lumM(f),x,'line','g-o');
loghist(lumB(f),x,'line','g--');
loghist(starM(f),x,'line','b-');
legend('DynMass','DynA','DynFree','DynBias','LumMass','LumML','LumBias','Star');
xlabel('Mass [Msun/h]');
ylabel('Count');
print('-depsc','MLWLMassDistributionCorrected.eps');
%%
f=grp.Mult>2;
powerlaw_err=@(x,lnAerr,Berr,corrAB) sqrt(lnAerr^2+(log(x)*Berr).^2+2*log(x)*corrAB*lnAerr*Berr);
powerlawmodel_err=@(m,perr,c) powerlaw_err(m/1e12,perr(1)/log10(exp(1)),perr(2),c); %log(error) for y=10^p(1)*(m/1e12).^p(2)
powerlawmodel=@(m,p) 10^p(1)*(m/1e12).^p(2);
p2err=[0.29,0.14];
cov2=-0.9;
p3err=[0.35,0.18];
cov3=-0.95;
figure;
x=logspace(12,15);
h1=semilogx(grp.LumMass(f),lumB(f)./grp.LumMass(f),'r.');
hold on;
y=powerlawmodel(x,p3);
yerrln=powerlawmodel_err(x,p3err,cov3);
plot(x,exp(log(y)+yerrln),'r--');
plot(x,exp(log(y)-yerrln),'r--');
hold on;
h2=semilogx(grp.MassProxy(f),dynB(f)./grp.MassProxy(f),'g.');
y=powerlawmodel(x,p2);
yerrln=powerlawmodel_err(x,p2err,cov2);
plot(x,exp(log(y)+yerrln),'g--');
plot(x,exp(log(y)-yerrln),'g--');
f10=grp.Mult>9;
p3N10err=[0.5,0.24];cov3N10=-0.98;
h3=semilogx(grp.LumMass(f10),lumB0N10(f10)./grp.LumMass(f10),'k.');
x=logspace(13,15);
y=powerlawmodel(x,p3N10);
yerrln=powerlawmodel_err(x,p3N10err,cov3N10);
plot(x,exp(log(y)+yerrln),'k--');
plot(x,exp(log(y)-yerrln),'k--');
ylim([0,2]);
xlim([1e10,1e16]);
xlabel('LumMass or DynMass');
ylabel('Mass Bias');
h4=plot(xmed,ymed,'bo-');
legend([h1,h2,h3,h4],'LumMass','DynMass','LumMass (N>9)','MockDyn');
% h1=semilogx(grp.MassProxy(f),lumB(f)./grp.LumMass(f),'rx');

print('-depsc','MassBias.eps');
%%
figure;
nbin=10;
f=grp.Mult>2;
l=grp.TotFluxProxy(f);
[xm,ym]=skeleton(l,dynB(f)./l,nbin,0.683);
loglog(xm,ym,'r');
hold on;
[xm,ym,ybrd]=skeleton(l,lumM(f)./l,nbin,0.683);
loglog(xm,ym,'g');
[xm,ym,ybrd]=skeleton(l,starM(f)./l,nbin,0.683);
loglog(xm,ym,'k');
% ploterr(xm,ym,[],{ylim(:,1),ylim(:,2)},'r','logxy');
%%
myfigure;
ml=[9.2,2.265;9.5,2.2;9.8,2.08;10.1,2.035;10.4,1.95;10.7,2.07;11,2.22;11.3,2.4;11.6,2.48;11.88,2.55];
ml=10.^ml;
loglog(ml(:,1),ml(:,2),'s-');hold on;
% loglog(l,dynA(f)./l,'r.');
% hold on;
% loglog(l,dynF(f)./l,'g.');hold on;
% plot(l,dynB(f)./l,'b.');
% plot(l,lumB(f)./l,'c.');
loglog(l,lumM(f)./l,'m.');hold on;
% plot(l,starM(f)./l,'k.');
% loglog(l,grp.LumMass(f)./l,'r.');hold on;
ylim([1,1e3])
A=1.;
% A=10^0.2;
ploterr([10^10.5,2e11,5e11,10^12.5],[80,132,177,83]*A,{[1e10,1e11,3e11,1e12],[1e11,3e11,1e12,1e13]},[42,45,57,77],'ko','logxy')
p4err=[0.1,0.22];
% x=logspace(10,13);
% yu=10^(p4(1)+p4err(1))*(x/1e12).^(p4(2)+p4err(2));
% yl=10^(p4(1)-p4err(1))*(x/1e12).^(p4(2)-p4err(2));
% plot(x,yu./x,'-',x,yl./x,'-');
xlabel('$L_r/(h^{\frac{\ }{\ }2}L_{r,{\odot}})$');ylabel('$M/L_r [h M_\odot/L_{r,\odot}]$');
% legend([h1(1),h2(1),h3],'WL','DynMass','2PIGG','location','northwest');
% print('-depsc','ML_MLWL.eps');
%%
%LumMass binned fit
m=[1e13,3e13; 3e13,1e14; 1e14,1e15; 1e15,1.5e15]; %lummass bin
mfitval=10.^[13.04,13.53,13.98,14.83]'; %MLWL fitted val
mfit=10.^[13.04-0.18,13.04+0.18;13.53-0.13,13.53+0.13; 13.98-0.16,13.98+0.16; 14.83-0.41,14.83+0.41]; %error bounds
% patch(mm,mmfit)
f=grp.Mult>2;
mval=m(:,1);
merr=mval;
for i=1:size(m,1)
    mval(i)=median(grp.LumMass(f&grp.LumMass>m(i,1)&grp.LumMass<m(i,2)));
%     mval(i)=mean(grp.LumMass(f&grp.LumMass>m(i,1)&grp.LumMass<m(i,2)));
    merr(i)=std(grp.LumMass(f&grp.LumMass>m(i,1)&grp.LumMass<m(i,2)));
end
figure();
loglog(grp.LumMass(f),lumM0(f),'g.','markersize',5);hold on;
loglog(grp.LumMass(f),lumB0(f),'b.');
% loglog(grp.LumMass(f),lumB(f),'b.','markersize',1);
ploterr(mval,mfitval,{m(:,1),m(:,2)},{mfit(:,1),mfit(:,2)},'ro','logxy');
% ploterr(mval,mfitval,merr,{mfit(:,1),mfit(:,2)},'ro','logxy');
% mm=m(:,[1,2,2,1])';
% mmfit=mfit(:,[1,1,2,2])';
% plot_boxes(m,mfit,'r')
x=logspace(12,15);
p3err=[0.35,0.18];
y=powerlawmodel(x,p3).*x;
yerrln=powerlawmodel_err(x,p3err,cov3);
plot(x,exp(log(y)+yerrln),'b--');
plot(x,exp(log(y)-yerrln),'b--');
plot([1e11,1e16],[1e11,1e16],'k--')
xlabel('LumMass');
ylabel('LumBias Mass (Raw fit)');
print('-depsc','LumMassFitRaw.eps')