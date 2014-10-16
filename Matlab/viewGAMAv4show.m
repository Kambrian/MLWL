%%
% view GAMA catalogue G3Cv4
%%
cd /work/Projects/Lensing/data
grpmock=loadGAMAcsv('G3Cv4/mock/G3CMock1FoFGroup194v04.dat',4);
halomock=loadGAMAcsv('G3Cv4/mock/G3CMock1HaloGroup194v04.dat',4);
% halomock.Mult=halomock.Nhalo;halomock=rmfield(halomock,'Nhalo');
% halomock.Zfof=halomock.Zhalo;halomock=rmfield(halomock,'Zhalo');
galmock=loadGAMAcsv('G3Cv4/mock/G3CMock1Galv04.dat',4);
% galmockv6=loadGAMAcsv('201309/DMUG3Cv06/mocks/G3CMockGalv06.csv',6);
refmock=loadGAMAcsv('G3Cv4/mock/G3CMock1Ref194v04.dat',4);
Ac=-1.2;An=20.7;Az=2.3;
Bc=0.94;Bn=-0.67;Bz=0.16;
halomock.MassProxy=halomock.MassProxy.*(Ac+An./sqrt(halomock.Mult)+Az./sqrt(halomock.Zfof));
halomock.TotFluxProxy=halomock.TotFluxProxy.*(Bc+Bn./sqrt(halomock.Mult)+Bz./sqrt(halomock.Zfof));
[halomock.LumMass,C1]=luminosity_mass(halomock);
grpmock.MassProxyRaw=grpmock.MassProxy;
grpmock.TotFluxProxyRaw=grpmock.TotFluxProxy;
grpmock.MassProxy=grpmock.MassProxy.*(Ac+An./sqrt(grpmock.Mult)+Az./sqrt(grpmock.Zfof));
grpmock.TotFluxProxy=grpmock.TotFluxProxy.*(Bc+Bn./sqrt(grpmock.Mult)+Bz./sqrt(grpmock.Zfof));
[grpmock.LumMass,C2]=luminosity_mass(grpmock);

grp=loadGAMAcsv('G3Cv4/group/G3CFoFGroup194v04.dat',4);
gal=loadGAMAcsv('G3Cv4/group/G3CGalv04.dat',4);
ref=loadGAMAcsv('G3Cv4/group/G3CRef194v04.dat',4);
grp.MassProxyRaw=grp.MassProxy;
grp.TotFluxProxyRaw=grp.TotFluxProxy;
grp.MassProxy=grp.MassProxy.*(Ac+An./sqrt(grp.Mult)+Az./sqrt(grp.Zfof));
grp.TotFluxProxy=grp.TotFluxProxy.*(Bc+Bn./sqrt(grp.Mult)+Bz./sqrt(grp.Zfof));
[grp.LumMass,C3]=luminosity_mass(grp);
[grp09,grp12,grp15]=split_mockcat(grp);

linkmock=[refmock.GalID,refmock.GroupID];  %galID: index (row number) for galaxies in the galaxy file
linkmock=sortrows(linkmock,2);
link=[ref.GalID,ref.GroupID];
link=sortrows(link,2);

cd /work/Projects/Lensing/outputv4
%%
figure;plot([C1(:),C2(:),C3(:)]);
figure;loglog(halomock.HaloMass,halomock.MassProxy,'.');
hold on;
loglog(halomock.HaloMass,halomock.LumMass,'r.');
plot([1e10,1e15],[1e10,1e15]);
% loglog(halomock.HaloMass,halomock.HybMass,'g.');
%%
figure;loglog(grpmock.LumMass,grpmock.MassProxy,'.');
hold on;
plot([1e10,1e15],[1e10,1e15]);

figure;loglog(grp.LumMass,grp.MassProxy,'.');
hold on;
plot([1e10,1e15],[1e10,1e15]);

figure;loglog(grp.LumMass,grp.TotFluxProxy,'.');
figure;semilogx(grp.TotFluxProxy,grp.LumMass./grp.TotFluxProxy,'.')

%%
galmock.HaloMass=galmock.HaloID;
for i=1:numel(galmock.HaloID)
    j=find(halomock.HaloID==galmock.HaloID(i));
    if j
        galmock.HaloMass(i)=halomock.HaloMass(j);
    else
        galmock.HaloMass(i)=1e9;  %this should be 0, but use 1e9 to be visible in logplot
    end
end
%% galmock halo mass
mockmass.Complexity=grpmock.GroupID;  %number of different group masses
mockmass.Mmean=grpmock.GroupID;   % average mass
mockmass.Gmean=grpmock.GroupID;   % geometric average mass, or average in logspace
mockmass.Mmode=grpmock.GroupID;   % most frequent mass
mockmass.Purity=grpmock.GroupID;  % frequency of Mmode
mockmass.MBCG=galmock.HaloMass(grpmock.BCGRef); %group mass for BCG
mockmass.MIter=galmock.HaloMass(grpmock.IterCenRef); %group mass for BCG
GrpOffset=sum(~linkmock(:,2));
for gid=1:numel(grpmock.GroupID)
    mass=galmock.HaloMass(linkmock(GrpOffset+(1:grpmock.Mult(gid)),1));
    ids=galmock.HaloID(linkmock(GrpOffset+(1:grpmock.Mult(gid)),1));
    ids=[ids,galmock.GalID(linkmock(GrpOffset+(1:grpmock.Mult(gid)),1))];
    GrpOffset=GrpOffset+grpmock.Mult(gid);
    ids=sortrows(ids);
    mockmass.Complexity(gid)=sum(logical(diff(ids(:,1))))+1;
    mockmass.Mmean(gid)=mean(mass);
    mockmass.Gmean(gid)=geomean(mass(logical(mass)));  %geomean of non-zero elements
    [mockmass.Mmode(gid),mockmass.Purity(gid)]=mode(ids(:,1));
    mid=ids(ids(:,1)==mockmass.Mmode(gid),2);mid=mid(1); %take GalID of the mode galaxy
    mockmass.Mmode(gid)=galmock.HaloMass(mid);
    if mockmass.Purity(gid)==1
        mockmass.Mmode(gid)=max(mass);
    end
end
mockmass.Purity=mockmass.Purity./grpmock.Mult;
%%
nq=max(grpmock.Mult);
Q=zeros(nq,1);eQ=Q;
for i=1:nq
    Q(i)=mean(mockmass.Purity(grpmock.Mult==i));
    eQ(i)=std(mockmass.Purity(grpmock.Mult==i));
end
x=-0.05:0.1:1.05;
myfigure;
[xm,ym]=linhist(mockmass.Purity,x);
f2=grpmock.Mult==2;
f3=grpmock.Mult==3;
f4=grpmock.Mult>=4;
[xm2,ym2]=linhist(mockmass.Purity(f2),x);
[xm3,ym3]=linhist(mockmass.Purity(f3),x);
[xm4,ym4]=linhist(mockmass.Purity(f4),x);
bar(0:0.1:1,[ym;ym2;ym3;ym4]','hist');
set(gca,'xtick',0:0.1:1);
xlim([0,1.05]);
xlabel('Mock Group Purity');ylabel('Counts');
l=legend('All','Mult=2','Mult=3','Mult>3');set(l,'location','northwest');
print('-depsc','MockGroupPurity.eps');
myfigure;errorbar(1:nq,Q,eQ,'.');
ylabel('Average Purity');xlabel('Multiplicity');
xlim([0,nq+1]);ylim([0.4,1.1]);
set(gca,'xscale','log');
% print('-depsc','MockGroupQuality.eps');
%% Complexity Histograms
x=min(mockmass.Complexity)-0.5:max(mockmass.Complexity)+0.5;
myfigure;[xm,ym]=linhist(mockmass.Complexity,x);
f2=grpmock.Mult==2;
f3=grpmock.Mult==3;
f4=grpmock.Mult>=4;
[xm2,ym2]=linhist(mockmass.Complexity(f2),x);
[xm3,ym3]=linhist(mockmass.Complexity(f3),x);
[xm4,ym4]=linhist(mockmass.Complexity(f4),x);
bar(xm,[ym;ym2;ym3;ym4]');
set(gca,'yscale','log');
set(gca,'xtick',0:15)
xlabel('Mock Group Complexity');ylabel('Counts');
legend('All','Mult=2','Mult=3','Mult>3');
print('-depsc','MockMassComplexity.eps');
%%
figure;loglog(mockmass.MBCG,mockmass.Mmode,'.');
figure;loglog(mockmass.MBCG,mockmass.Mmean,'.');
f=grpmock.Mult==2;

%%
f=grpmock.Mult<=2;
myfigure;loglog(mockmass.Mmode(f),mockmass.Mmean(f),'g.','markersize',5);
f=grpmock.Mult>2;
hold on;loglog(mockmass.Mmode(f),mockmass.Mmean(f),'ro');
legend('Mult<=2','Mult>2');legend('location','southeast');
xlabel('Mass Mode');ylabel('Mass Mean');
print('-depsc','MockMassCmpMean_Mode.eps');
%% Compare BCG mass and Mmode
fBCG=mockmass.MBCG==mockmass.Mmode;
x=1:max(grpmock.Mult);
y=zeros(size(x));
for i=x
    y(i)=sum(fBCG(grpmock.Mult==i))/numel(fBCG(grpmock.Mult==i));
end
myfigure;plot(x,y,'o');
xlabel('Multiplicity');ylabel('BCG-Meaningful Fraction');
hold on;
plot(x,repmat(sum(fBCG)/numel(fBCG),size(x)),'-');
plot(x,Q,'r.');
set(gca,'xscale','log');
l=legend('BCG','Total','Purity');
set(l,'location','southeast');
% print('-depsc','MockGroupBCGRepFraction.eps');
% figure;loglog(mockmass.MBCG(~fBCG),mockmass.Mmode(~fBCG),'.')
%% Compare Iter mass and Mmode
fIter=mockmass.MIter==mockmass.Mmode;
x=1:max(grpmock.Mult);
y=zeros(size(x));
for i=x
    y(i)=sum(fIter(grpmock.Mult==i))/numel(fIter(grpmock.Mult==i));
end
myfigure;plot(x,y,'o');
xlabel('Multiplicity');ylabel('Iter-Meaningful Fraction');
hold on;
plot(x,repmat(sum(fIter)/numel(fIter),size(x)),'-');
plot(x,Q,'r.');
set(gca,'xscale','log');
l=legend('Iter','Total','Purity');
set(l,'location','southeast');
% print('-depsc','MockGroupIterRepFraction.eps');
% figure;loglog(mockmass.MIter(~fIter),mockmass.Mmode(~fIter),'.')
%% Purity histogram
f=grpmock.Mult>=2;
myfigure;linhist(mockmass.Purity(f),20,'stairsnorm');xlabel('galmock Group Mass Purity');ylabel('Fraction');
% print('-depsc','MockMassPurity.eps');
%% Mass distributions
myfigure;loghist(mockmass.Mmean,30,'stairs','r');
hold on;loghist(mockmass.Mmode,30,'stairs','g');
loghist(mockmass.MBCG,30,'stairs','b');
loghist(mockmass.MIter,30,'stairs','k');
legend('Mean','Mode','BCG','Iter');xlabel('Mass');ylabel('Counts');
print('-depsc','MockMassDistrMean_Mode.eps');
%% Min and Max mass for 2-member
f=grpmock.Mult<=2;
myfigure;loglog(mockmass.Mmode(f),2*mockmass.Mmean(f)-mockmass.Mmode(f),'b.','markersize',5);
xlabel('Mmax');ylabel('Mmin');
title('Mult=2');
print('-depsc','MockMassMinMax.eps')
%% mass observables
f=grpmock.Mult>2;
% figure;loglog(mockmass.Mmode(f),grpmock.MassProxy(f),'.');
% figure;loglog(mockmass.Mmean(f),grpmock.MassProxy(f),'.');
myfigure;
loglog(mockmass.MBCG(f),grpmock.LumMass(f),'r.','markersize',5);
hold on;
loglog(mockmass.MBCG(f),grpmock.MassProxy(f),'g.','markersize',5);hold on;
% loglog(mockmass.MBCG(f),sqrt(grpmock.MassProxy(f).*grpmock.LumMass(f)),'b.','markersize',5);
% xlim([1e10,1e17]);ylim([1e10,1e17]);
legend('Luminosity Mass','Dynamical Mass','Hybrid Mass');legend('location','northwest');
xlabel('$M_{BCGhost}/(M_{\odot}/h)$');ylabel('$Mass(M_{\odot}/h)$');title('$Mult>2$');
plot([1e10,1e17],[1e10,1e17],'r-');
% print('-depsc','MassObservable.eps');
%%
mtrue=mockmass.MIter;
% mtrue=max([mockmass.MIter,mockmass.Mmode],[],2);
mobs=grpmock.MassProxy;
f=grpmock.Mult>2&mtrue>1e1;%&mockmass.Purity>0.0&mockmass.Complexity<10;
x=logspace(10,16,10);
figure();
% semilogx(mobs(f),mtrue(f)./mobs(f),'.');hold on;
[xmed,ymed,yl]=skeleton(mobs(f),mtrue(f)./mobs(f),x,0.68);
semilogx(xmed,ymed,'rd-');hold on;
plot(xm,ym,'g');hold on;
plot(xmed,yl(:,1),'r--');
plot(xmed,yl(:,2),'r--');
ylim([0,1.5])
%%
mtrue=mockmass.MIter;
% mtrue=min(mockmass.MIter,mockmass.Gmean);
mobs=grpmock.LumMass;
f=grpmock.Mult>2&mtrue>1e1;%&mockmass.Purity>0.0&mockmass.Complexity<10;
x=logspace(10,16,20);
figure();
% semilogx(mobs(f),mtrue(f)./mobs(f),'.');hold on;
[xmed,ymed,yl]=skeleton(mobs(f),mtrue(f)./mobs(f),x,0.68);
semilogx(xmed,ymed,'rd-');hold on;
plot(xm,ym,'g');hold on;
plot(xmed,yl(:,1),'r--');
plot(xmed,yl(:,2),'r--');
ylim([0,1.5])
%%
figure;
f=grpmock.Mult>2;
loglog(grpmock.VelDisp(f).^2,mtrue(f),'.')
x=logspace(log10(100),log10(2000),6);
[xmed,ymed,yl,xm,ym,yerr]=skeleton(grpmock.VelDisp(f),log10(mockmass.MIter(f)),x,0.68);
ploterr(xmed,ymed,[],yerr,'rd-','logx');hold on;
%%
nbin=30;
f=grpmock.Mult>2;
myfigure;
[xx,yy,n,s]=densitygrid(log10(mockmass.MIter(f)),log10(grpmock.LumMass(f)),[nbin,nbin],[10,17],[10,17]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(grpmock.VelDisp'),minmax(grpmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
hold on;plot([11,16],[11,16],'w');
[xmed,ymed]=skeleton(log10(grpmock.LumMass(f)),log10(mockmass.MIter(f)),10,0.68);
plot(ymed,xmed,'w');
xlabel('$IterMass$');
ylabel('$LumMass$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
nbin=30;
myfigure;
[xx,yy,n,s]=densitygrid(log10(mockmass.MIter(f)),log10(grpmock.MassProxy(f)),[nbin,nbin],[10,17],[10,17]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(grpmock.VelDisp'),minmax(grpmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
hold on;plot([11,16],[11,16],'w');
[xmed,ymed]=skeleton(log10(grpmock.MassProxy(f)),log10(mockmass.MIter(f)),10,0.68);
plot(ymed,xmed,'w');
xlabel('$IterMass$');
ylabel('$DynMass$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
nbin=30;
myfigure;
[xx,yy,n,s]=densitygrid(log10(mockmass.MIter(f)),log10(sqrt(grpmock.MassProxy(f).*grpmock.LumMass(f))),[nbin,nbin],[10,17],[10,17]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(grpmock.VelDisp'),minmax(grpmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
hold on;plot([11,16],[11,16],'w');
[xmed,ymed]=skeleton(log10(sqrt(grpmock.MassProxy(f).*grpmock.LumMass(f))),log10(mockmass.MIter(f)),10,0.68);
plot(ymed,xmed,'w');
xlabel('$IterMass$');
ylabel('$HybMass$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
%%
f=grpmock.Mult==2;
% figure;loglog(mockmass.Mmode(f),grpmock.MassProxy(f),'.');
% figure;loglog(mockmass.Mmean(f),grpmock.MassProxy(f),'.');
myfigure;loglog(mockmass.MBCG(f),grpmock.LumMass(f),'r.','markersize',5);hold on;loglog(mockmass.MBCG(f),grpmock.MassProxy(f),'g.','markersize',5);
% loglog(mockmass.MBCG(f),sqrt(grpmock.MassProxy(f).*grpmock.LumMass(f)),'b.','markersize',5);
% xlim([1e10,1e17]);ylim([1e10,1e17]);
legend('Luminosity Mass','Dynamical Mass','Hybrid Mass');legend('location','northwest');
xlabel('$M_{BCGhost}/(M_{\odot}/h)$');ylabel('$Mass(M_{\odot}/h)$');title('$Mult=2$');
plot([1e10,1e17],[1e10,1e17],'r-');
% print('-depsc','MassObservable2.eps');
%% mass observable--best construction
mockmass.Mbest=grpmock.GroupID;
mockmass.Mbest(grpmock.Mult<3)=grpmock.LumMass(grpmock.Mult<3);
mockmass.Mbest(grpmock.Mult>=3)=sqrt(grpmock.LumMass(grpmock.Mult>=3).*grpmock.MassProxy(grpmock.Mult>=3));
myfigure;loglog(mockmass.MBCG(f),mockmass.Mbest(f),'b.','markersize',4);hold on;
xlim([1e10,1e17]);ylim([1e10,1e17]);
xlabel('$M_{BCGhost}/(M_{\odot}/h)$');ylabel('$Mbest(M_{\odot}/h)$');
title({'$M_{best}=M_{Lum}~(Mult=2),=\sqrt(M_{Lum}M_{dyn})~(Mult>=3)$'})
plot([1e10,1e17],[1e10,1e17],'r-');
% print('-depsc','MassObservable3.eps');
%%
[xx,yy,n,s]=densitygrid(mockmass.Complexity,mockmass.Purity,[16,20]);
figure;pcolor(xx,yy,log10(n+1));xlabel('Complexity');ylabel('Purity');colormap('gray');
[xx,yy,n,s]=densitygrid(grpmock.Mult,mockmass.Purity,[153,20]);
figure;pcolor(xx,yy,log10(n+1));xlabel('Multiplicity');ylabel('Purity');colormap('gray');
[xx,yy,n]=gridcounts(grpmock.Mult,mockmass.Complexity);
figure;pcolor(xx-0.5,yy-0.5,log10(n+1));xlabel('Multiplicity');ylabel('Complexity');colormap('gray');
% print('-depsc','MockGroupComplexity-MultJoint.eps');
%% LumMass is the best mass
f=grpmock.Mult>2;
% figure;loglog(mockmass.Mmode(f),grpmock.MassProxy(f),'.');
% figure;loglog(mockmass.Mmean(f),grpmock.MassProxy(f),'.');
myfigure;loglog(mockmass.Mmode(f),grpmock.MassProxy(f),'r.','markersize',5);hold on;loglog(mockmass.Mmode(f),grpmock.LumMass(f),'g.','markersize',5);
% loglog(mockmass.Mmode(f),sqrt(grpmock.MassProxy(f).*grpmock.LumMass(f)),'b.','markersize',5);
% xlim([1e10,1e17]);ylim([1e10,1e17]);
legend('Dynamical Mass','Luminosity Mass','Combined Mass');legend('location','northwest');
xlabel('$M_{halo}/(M_{\odot}/h)$');ylabel('$Mass(M_{\odot}/h)$');title('$Mult>2$');
plot([1e10,1e17],[1e10,1e17],'r-');
% print('-depsc','MassObservable.eps');
myfigure;loglog(mockmass.Mmode(f),mockmass.Mbest(f),'b.','markersize',5);hold on;
loglog(mockmass.Mmode(f),grpmock.LumMass(f),'k.','markersize',5);
xlim([1e10,1e17]);ylim([1e10,1e17]);
xlabel('$M_{halo}/(M_{\odot}/h)$');ylabel('$Mass(M_{\odot}/h)$');title('$Mult>2$');
plot([1e10,1e17],[1e10,1e17],'r-');
% figure;loglog(grpmock.LumMass(f),grpmock.MassProxy(f),'.');
% figure;loglog(mockmass.Mmode(f),(grpmock.MassProxy(f)+grpmock.LumMass(f))/2,'b.','markersize',5);xlim([1e10,1e17]);ylim([1e10,1e17]);
% figure;loglog(mockmass.Mmean(f),grpmock.TotFluxInt(f),'.');
% figure;loglog(mockmass.Mmean(f),grpmock.TotFlux(f),'.');
% figure;loglog(mockmass.Mmean(f),grpmock.BCGmag(f),'.');
% figure;loglog(mockmass.Mmean(f),grpmock.VelDisp(f),'.');
% figure;loglog(mockmass.Mmean(f),grpmock.Rad1Sig(f),'.');
%%
f=grpmock.Mult>2;
[xm,ym]=skeleton(log10(mockmass.Mmode(f)),log10(grpmock.LumMass(f)),10,0.68);
[xm2,ym2]=skeleton(log10(grpmock.LumMass(f)),log10(mockmass.Mmode(f)),10,0.68);
plot(log10(mockmass.Mmode(f)),log10(grpmock.LumMass(f)),'.');
hold on;
plot(xm,ym,'ro-');
plot(ym2,xm2,'gd-');
plot([11,16],[11,16],'k-')
%% Mass estimator dispersion:
f=grpmock.Mult>2;
% figure;
% loglog(grpmock.MassProxyRaw(f),mockmass.Mmode(f),'r.');
figure;
loglog(grpmock.TotFluxProxy(f),mockmass.Mmode(f),'g.');
%%
m=logspace(12,15,6);
mu=zeros(5,1);
sigma=mu;
mumu=[];
for i=1:5
%     mo=mockmass.Mmode;
    mo=grpmock.LumMass;
%     mo=grpmock.MassProxy;
    f=mo>m(i)&mo<m(i+1);
    f=f&grpmock.Mult>2;
    f=f&mockmass.Mmode>1e10;
    f=f&grpmock.MassProxy>1e10;
%     mu(i)=log10(median(mo(f)));
%     [mu(i),sigma(i)]=normfit(log10(mockmass.Mmode(f)));
%     [mu(i),sigma(i)]=normfit(log10(grpmock.MassProxy(f)))
%     mu(i)=median(log10(grpmock.LumMass(f)));
    mu(i)=median(log10(mockmass.Mmode(f)));
    mumu=[mumu;log10(mean(mo(f))),median(log10(mo(f)))];
end
figure;
plot(mumu(:,1),mu,'r');
hold on;
plot(mumu(:,2),mu,'g');
plot([12,15],[12,15],'-');

%%
myfigure;
plot(mu,sigma,'d-')
xlabel('log10(M)');ylabel('$\sigma (log10(M))$');
title('Dynamical Mass Bin');
% print('-depsc','MassScatter_DynMass.eps');
%% group redshift: not worth examination;
mockzz.Mean=grpmock.GroupID;   % average mass
mockzz.Med=grpmock.GroupID;   % most frequent mass
mockzz.Disp=grpmock.GroupID;  % frequency of Mmode
mockzz.BCG=grpmock.GroupID;
mockzz.Ext=grpmock.GroupID;
GrpOffset=sum(~linkmock(:,2));
for gid=1:numel(grpmock.GroupID)
    zz=galmock.Zspec(linkmock(GrpOffset+(1:grpmock.Mult(gid)),1));
    mockzz.BCG(gid)=zz(linkmock(GrpOffset+(1:grpmock.Mult(gid)),1)==grpmock.BCGRef(gid));
    mockzz.Mean(gid)=mean(zz);
    mockzz.Med(gid)=median(zz);
    mockzz.Disp(gid)=std(zz);
    mockzz.Ext(gid)=max(zz)-min(zz);
    GrpOffset=GrpOffset+grpmock.Mult(gid);
end
figure;
% plot(mockzz.Med,mockzz.Mean,'r.');
% hold on;
plot(mockzz.Med,mockzz.BCG,'g.','markersize',2);
figure;
plot(mockzz.Med,mockzz.Disp,'bo');hold on;
plot(mockzz.Med,mockzz.Ext,'ko');
figure;
hist(mockzz.Med-mockzz.BCG,20)
%% haloMass+Redshift joint distribution
nbin=30;
% f=mockmass.Mmode>0;
f=ones(size(grpmock.Mult));f=logical(f);
% f=f&grpmock.Mult>=2;
mass=grpmock.LumMass;
% mass=grpmock.MassProxy;
myfigure;
[xx,yy,n,s]=densitygrid(log10(mass(f)),grpmock.Zfof(f),[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$log(M_{Lum}/(M_{\odot}/h))$');
ylabel('$z$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','mockGroupMassRedshift.eps');
%%
nbin=5;
%% Multiplicity fun
myfigure;
% [xm,ym,dym]=loghist(grpmock.Mult,nbin,'stairs','k--');
% hold on;
% [x1,y1,dy1]=loghist(cat09.Mult,nbin,'stairs','r');
% [x2,y2,dy2]=loghist(cat12.Mult,nbin,'stairs','g');
% [x3,y3,dy3]=loghist(cat15.Mult,nbin,'stairs','b');
[x,y,dy]=loghist(grp.Mult(grp.Mult>1),nbin,'stairs','k');
xlabel('Group Multiplicity');
ylabel('dN');
hold on;
x=[3,4,5,6,7,8,9,12,18,26,41,71,220];
y=[58788,27083,14925,8744,5630,3858,6196,4427,1711,787,272,47,0];
hold on;stairs(x,y,'m');
% legend('galmock','G09','G12','G15','GAMA','Jhonston','location','best');
legend('GAMA','Jhonston','location','best');
print('-depsc','GroupMultfunCmp.eps');

%% Massfun
myfigure;
[xm,ym,dym]=loghist(grpmock.MassProxy,nbin,'stairs','k--');
hold on;
[x1,y1,dy1]=loghist(cat09.Mass,nbin,'stairs','r');
[x2,y2,dy2]=loghist(cat12.Mass,nbin,'stairs','g');
[x3,y3,dy3]=loghist(cat15.Mass,nbin,'stairs','b');
[x,y,dy]=loghist([cat09.Mass;cat12.Mass;cat15.Mass],nbin,'stairs','k');
xlabel('$M_{dyn}/(M_{\odot}/h)$');
ylabel('dN');
legend('galmock','G09','G12','G15','GAMA','location','best');
print('-depsc','GroupMassfun.eps');

dym=dym/sum(ym);dy1=dy1./sum(y1);dy2=dy2./sum(y2);dy3=dy3/sum(y3);dy=dy/sum(y);
myfigure;
loglog(xm,dym,'k.--');
hold on;
loglog(x1,dy1,'r.-');
loglog(x2,dy2,'g.-');
loglog(x3,dy3,'b.-');
loglog(x,dy,'k.-');
xlabel('$M_{dyn}/(M_{\odot}/h)$');
ylabel('$dN/dM_{dyn}/N_{tot}$');
legend('galmock','G09','G12','G15','GAMA','location','best');

print('-depsc','GroupMassfunSpec.eps');
%% Stellar Mass Fun
myfigure;
[x1,y1,dy1]=loghist(cat09.TotSMInt,nbin,'stairs','r');
hold on;
[x2,y2,dy2]=loghist(cat12.TotSMInt,nbin,'stairs','g');
[x3,y3,dy3]=loghist(cat15.TotSMInt,nbin,'stairs','b');
[x,y,dy]=loghist(cat.TotSMInt,nbin,'stairs','k');
xlabel('$M_{star}/(M_{\odot}/h)$');
ylabel('dN');
legend('G09','G12','G15','GAMA','location','best');

print('-depsc','GroupStellarfun.eps');

dy1=dy1./sum(y1);dy2=dy2./sum(y2);dy3=dy3/sum(y3);dy=dy/sum(y);
myfigure;
loglog(x1,dy1,'r.-');
hold on;
loglog(x2,dy2,'g.-');
loglog(x3,dy3,'b.-');
loglog(x,dy,'k.-');
xlabel('$M_{star}/(M_{\odot}/h)$');
ylabel('$dN/dM_{star}/N_{tot}$');
legend('G09','G12','G15','GAMA','location','best');

print('-depsc','GroupStellarfunSpec.eps');
%% LumMassfun
myfigure;
[xm,ym,dym]=loghist(grpmock.LumMass,nbin,'stairs','k--');
hold on;
[x1,y1,dy1]=loghist(cat09.LumMass,nbin,'stairs','r');
[x2,y2,dy2]=loghist(cat12.LumMass,nbin,'stairs','g');
[x3,y3,dy3]=loghist(cat15.LumMass,nbin,'stairs','b');
[x,y,dy]=loghist(cat.LumMass,nbin,'stairs','k');
xlabel('$M_{lum}/(M_{\odot}/h)$');
ylabel('dN');
legend('galmock','G09','G12','G15','GAMA','location','best');
print('-depsc','GroupLumMassfun.eps');

dym=dym/sum(ym);dy1=dy1./sum(y1);dy2=dy2./sum(y2);dy3=dy3/sum(y3);dy=dy/sum(y);
myfigure;
loglog(xm,dym,'k.--');
hold on;
loglog(x1,dy1,'r.-');
loglog(x2,dy2,'g.-');
loglog(x3,dy3,'b.-');
loglog(x,dy,'k.-');
xlabel('$M_{lum}/(M_{\odot}/h)$');
ylabel('$dN/dM_{lum}/N_{tot}$');
legend('galmock','G09','G12','G15','GAMA','location','best');

print('-depsc','GroupLumMassfunSpec.eps');
%% Lumosity Fun
myfigure;
[xm,ym,dym]=linhist(grpmock.TotFluxInt,nbin,'stairs','k--');
hold on;
[x1,y1,dy1]=linhist(cat09.TotFluxInt,nbin,'stairs','r');
[x2,y2,dy2]=linhist(cat12.TotFluxInt,nbin,'stairs','g');
[x3,y3,dy3]=linhist(cat15.TotFluxInt,nbin,'stairs','b');
[x,y,dy]=linhist(cat.TotFluxInt,nbin,'stairs','k');
set(gca,'yscale','log');
xlabel('$M_{r}$');
ylabel('dN');
legend('galmock','G09','G12','G15','GAMA','location','south');
print('-depsc','GroupLumfun.eps');

dym=dym/sum(ym);dy1=dy1./sum(y1);dy2=dy2./sum(y2);dy3=dy3/sum(y3);dy=dy/sum(y);
myfigure;
loglog(xm,dym,'k.--');
hold on;
loglog(x1,dy1,'r.-');
loglog(x2,dy2,'g.-');
loglog(x3,dy3,'b.-');
loglog(x,dy,'k.-');
set(gca,'xscale','linear');
xlabel('$M_{r}$');
ylabel('$dN/dM_{r}/N_{tot}$');
legend('galmock','G09','G12','G15','GAMA','location','south');
print('-depsc','GroupLumfunSpec.eps');
%% Redshift Distr
myfigure;
linhist(grpmock.Zfof,nbin,'stairs','k--');
hold on;
linhist(cat09.Zfof,nbin,'stairs','r');
linhist(cat12.Zfof,nbin,'stairs','g');
linhist(cat15.Zfof,nbin,'stairs','b');
linhist([cat09.Zfof;cat12.Zfof;cat15.Zfof],nbin,'stairs','k');
set(gca,'yscale','log');
xlabel('$z$');
ylabel('dN');
legend('galmock','G09','G12','G15','GAMA','location','best');

print('-depsc','GroupRedshift.eps');
%% Redshift Completeness
xbin=(0.05:0.1:1.05)';
myfigure;
linhist(grpmock.ZComp,xbin,'stairsnorm','k--');
hold on;
linhist(cat09.ZComp,xbin,'stairsnorm','r');
linhist(cat12.ZComp,xbin,'stairsnorm','g');
linhist(cat15.ZComp,xbin,'stairsnorm','b');
linhist(cat.ZComp,xbin,'stairsnorm','k');
xlabel('Completeness');
ylabel('Fraction');
legend('galmock','G09','G12','G15','GAMA','location','best');

print('-depsc','GroupZComplete.eps');
%% VD Distr
myfigure;
[xm,ym,dym]=loghist(grpmock.VelDisp,nbin,'stairs','k--');
hold on;
[x1,y1,dy1]=loghist(cat09.VelDisp,nbin,'stairs','r');
[x2,y2,dy2]=loghist(cat12.VelDisp,nbin,'stairs','g');
[x3,y3,dy3]=loghist(cat15.VelDisp,nbin,'stairs','b');
[x,y,dy]=loghist(cat.VelDisp,nbin,'stairs','k');
xlabel('$\sigma/(km/s)$');
ylabel('dN');
legend('galmock','G09','G12','G15','GAMA','location','best');
print('-depsc','GroupVelDisp.eps');

dym=dym/sum(ym);dy1=dy1./sum(y1);dy2=dy2./sum(y2);dy3=dy3/sum(y3);dy=dy/sum(y);
myfigure;
loglog(xm,dym,'k.--');
hold on;
loglog(x1,dy1,'r.-');
loglog(x2,dy2,'g.-');
loglog(x3,dy3,'b.-');
loglog(x,dy,'k.-');
% set(gca,'xscale','linear');
xlabel('$\sigma/(km/s)$');
ylabel('$dN/d\sigma/N_{tot}$');
legend('galmock','G09','G12','G15','GAMA','location','south');
print('-depsc','GroupVelDispSpec.eps');
%% Size Distr
myfigure;
[xm,ym,dym]=loghist(grpmock.Rad50,nbin,'stairs','k--');
hold on;
[x1,y1,dy1]=loghist(cat09.Rad50,nbin,'stairs','r');
[x2,y2,dy2]=loghist(cat12.Rad50,nbin,'stairs','g');
[x3,y3,dy3]=loghist(cat15.Rad50,nbin,'stairs','b');
[x,y,dy]=loghist(cat.Rad50,nbin,'stairs','k');
xlabel('$R_{0.5}/(Mpc/h)$');
ylabel('dN');
legend('galmock','G09','G12','G15','GAMA','location','best');
print('-depsc','GroupSize.eps');

dym=dym/sum(ym);dy1=dy1./sum(y1);dy2=dy2./sum(y2);dy3=dy3/sum(y3);dy=dy/sum(y);
myfigure;
loglog(xm,dym,'k.--');
hold on;
loglog(x1,dy1,'r.-');
loglog(x2,dy2,'g.-');
loglog(x3,dy3,'b.-');
loglog(x,dy,'k.-');
set(gca,'xlim',[1e-4,10],'ylim',[1e-4,1e2]);
xlabel('$R_{0.5}/(Mpc/h)$');
ylabel('$dN/dR_{0.5}/N_{tot}$');
legend('galmock','G09','G12','G15','GAMA','location','southwest');
print('-depsc','GroupSizeSpec.eps');
%% Color
myfigure;
[x1,y1,dy1]=linhist(g09.AB_u-g09.AB_r,nbin,'stairs','r');
hold on;
[x2,y2,dy2]=linhist(g12.AB_u-g12.AB_r,nbin,'stairs','g');
[x3,y3,dy3]=linhist(g15.AB_u-g15.AB_r,nbin,'stairs','b');
[x,y,dy]=linhist(gal.AB_u-gal.AB_r,nbin,'stairs','k');
set(gca,'yscale','log');
xlabel('$M_{r}$');
ylabel('dN');
legend('G09','G12','G15','GAMA','location','northeast');
% print('-depsc','GAMAColorfun.eps');

dy1=dy1./sum(y1);dy2=dy2./sum(y2);dy3=dy3/sum(y3);dy=dy/sum(y);
myfigure;
semilogy(x1,dy1,'r.-');
hold on;
semilogy(x2,dy2,'g.-');
semilogy(x3,dy3,'b.-');
semilogy(x,dy,'k.-');
% set(gca,'xscale','linear');
xlabel('$U_{AB}-R_{AB}$');
ylabel('$dN/d(U-R)/N_{tot}$');
legend('G09','G12','G15','GAMA','location','northeast');
print('-depsc','GAMAColorfunSpec.eps');
%% VD+Size joint distribution
nbin=30;
myfigure;
[xx,yy,n,s]=densitygrid(grpmock.VelDisp,grpmock.Rad50,[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(grpmock.VelDisp'),minmax(grpmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$\sigma/(km/s)$');
ylabel('$R_{0.5}/(Mpc/h)$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','mockGroupSizeVD.eps');

myfigure;
[xx,yy,n,s]=densitygrid(cat.VelDisp,cat.Rad50,[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(grpmock.VelDisp'),minmax(grpmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$\sigma/(km/s)$');
ylabel('$R_{0.5}/(Mpc/h)$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','GAMAGroupSizeVD.eps');
%% Size+Redshift joint distr
nbin=30;
myfigure;
[xx,yy,n,s]=densitygrid(grpmock.Rad50,grpmock.Zfof,[nbin,nbin]);
% pcolor(xx,yy,log10(n+1));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(grpmock.VelDisp'),minmax(grpmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
ylabel('$z$');
xlabel('$R_{0.5}/(Mpc/h)$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','mockGroupSizeRedshift.eps');

myfigure;
[xx,yy,n,s]=densitygrid(cat.Rad50,cat.Zfof,[nbin,nbin]);
% pcolor(xx,yy,log10(n+1));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(grpmock.VelDisp'),minmax(grpmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
ylabel('$z$');
xlabel('$R_{0.5}/(Mpc/h)$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','GAMAGroupSizeRedshift.eps');
%% VD+Redshift joint distribution
nbin=30;
myfigure;
[xx,yy,n,s]=densitygrid(grpmock.VelDisp,grpmock.Zfof,[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(grpmock.VelDisp'),minmax(grpmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$\sigma/(km/s)$');
ylabel('$z$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','mockGroupVDRedshift.eps');

myfigure;
[xx,yy,n,s]=densitygrid(cat.VelDisp,cat.Zfof,[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(grpmock.VelDisp'),minmax(grpmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$\sigma/(km/s)$');
ylabel('$z$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','GAMAGroupVDRedshift.eps');
%% Mass+Redshift joint distribution
nbin=30;
f=grpmock.MassProxy>0;
f=f&grpmock.Mult>=2;
myfigure;
[xx,yy,n,s]=densitygrid(log10(grpmock.MassProxy(f)),grpmock.Zfof(f),[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(grpmock.VelDisp'),minmax(grpmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$log(M_{dyn}/(M_{\odot}/h))$');
ylabel('$z$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
% print('-depsc','mockGroupMassRedshift_3.eps');

f=cat.Mass>0;
f=f&grp.Mult>=2;
myfigure;
[xx,yy,n,s]=densitygrid(log10(cat.Mass(f)),cat.Zfof(f),[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(grpmock.VelDisp'),minmax(grpmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$log(M_{dyn}/(M_{\odot}/h))$');
ylabel('$z$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
% print('-depsc','GAMAGroupMassRedshift_3.eps');
%% LumMass+Redshift joint distribution
nbin=30;
f=grpmock.LumMass>0;
f=f&grpmock.Mult>=3;
myfigure;
[xx,yy,n,s]=densitygrid(log10(grpmock.LumMass(f)),grpmock.Zfof(f),[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(grpmock.VelDisp'),minmax(grpmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$log(M_{Lum}/(M_{\odot}/h))$');
ylabel('$z$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','mockGroupLumMassRedshift_3.eps');

f=cat.LumMass>0;
f=f&grp.Mult>=3;
myfigure;
[xx,yy,n,s]=densitygrid(log10(cat.LumMass(f)),cat.Zfof(f),[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(grpmock.VelDisp'),minmax(grpmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$log(M_{Lum}/(M_{\odot}/h))$');
ylabel('$z$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','GAMAGroupLumMassRedshift_3.eps');
%% Color-Mag joint distribution
nbin=30;
% myfigure;
% [xx,yy,n,s]=densitygrid(galmock.AB_u-galmock.AB_r,galmock.AB_r,[nbin,nbin]);
% % pcolor(xx,yy,log10(1+n));hold on;
% contourf(xx,yy,log10((n+1)))
% % imagesc(minmax(grpmock.VelDisp'),minmax(grpmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
% colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
% xlabel('$U_{AB}-R_{AB}$');
% ylabel('$R_{AB}$');
% title([num2str(nbin),'$\times$',num2str(nbin)]);
% print('-depsc','mockColorMag.eps');

myfigure;
[xx,yy,n,s]=densitygrid(gal.AB_u-gal.AB_r,gal.AB_r,[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(grpmock.VelDisp'),minmax(grpmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$U_{AB}-R_{AB}$');
ylabel('$R_{AB}$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','GAMAColorMag.eps');
%% Mass-Mult
f=grpmock.MassProxy>0;
myfigure;
loglog(grpmock.MassProxy(f),grpmock.Mult(f),'.','markersize',5);
xlabel('$M_{dyn}/(M_{\odot}/h)$');
ylabel('$Multiplicity$');
print('-depsc','mockGroupMassMult.eps');

f=cat.Mass>0;
myfigure;
loglog(cat.Mass(f),grp.Mult(f),'.','markersize',5);
xlabel('$M_{dyn}/(M_{\odot}/h)$');
ylabel('$Multiplicity$');
print('-depsc','GAMAGroupMassMult.eps');
%%
myfigure;
loglog(grpmock.VelDisp,grpmock.MassProxy,'.','markersize',3);
xlabel('$\sigma/(km/s)$');
ylabel('$M_{dyn}/(M_{\odot}/h)$');
print('-depsc','mockVDMass.eps');

%%
myfigure;
subplot(3,1,1);
plot(g09.RA,g09.DEC,'.','markersize',5);
axis equal;
subplot(3,1,2);
plot(g12.RA,g12.DEC,'.','markersize',5);
axis equal;
subplot(3,1,3);
plot(g15.RA,g15.DEC,'.','markersize',5);
axis equal;
print('-depsc','GAMAholes.eps');
myfigure;
plot(cat.CenRA,cat.CenDEC,'.')
%% Mult-Redshift 
colors=repmat('rgbcmk',1,10);
nbin=[2,3,5,10,20,100];
zbin=[0,0.1,0.2,0.3,0.4,0.5];
for i=1:numel(nbin)-1
    linhist(grp.Zfof(grp.Mult>=nbin(i)&grp.Mult<nbin(i+1)),6,'norm',colors(i));hold on;
%     linhist(grp.Mult(grp.Zfof>=zbin(i)&grp.Zfof<zbin(i+1)),2:2:20,'norm',colors(i));hold on;
end
