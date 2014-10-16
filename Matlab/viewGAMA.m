%%
%view G3Cv2
%%
cd /home/kam/Projects/Lensing/data
% mock=loadGAMAcsv('mocksv2/mockcutcatvol1.csv');
catmock=loadGAMAcsv('G3Cv2/mock/mockgroupcatv2opt194lim194vol1.csv');
refmock=loadGAMAcsv('G3Cv2/mock/mockgrefsv2opt194lim194vol1.csv');
mock=loadGAMAcsv('G3Cv2/mock/mockcutcat194vol/mockcutcat194vol1.csv');

g09=loadGAMAcsv('G3Cv2/group/cutcatG09.csv');
cat09=loadGAMAcsv('G3Cv2/group/groupcatv2G09opt194lim194.csv');
ref09=loadGAMAcsv('G3Cv2/group/grefsv2G09opt194lim194.csv');

g12=loadGAMAcsv('G3Cv2/group/cutcatG12.csv');
cat12=loadGAMAcsv('G3Cv2/group/groupcatv2G12opt194lim194.csv');
ref12=loadGAMAcsv('G3Cv2/group/grefsv2G12opt194lim194.csv');

g15=loadGAMAcsv('G3Cv2/group/cutcatG15.csv');
cat15=loadGAMAcsv('G3Cv2/group/groupcatv2G15opt194lim194.csv');
ref15=loadGAMAcsv('G3Cv2/group/grefsv2G15opt194lim194.csv');

catv2=structcat([cat09;cat12;cat15]);
galv2=structcat([g09;g12;g15]);

linkmock=[refmock.galID,refmock.groupID];  %galID: index (row number) for galaxies in the galaxy file
linkmock=sortrows(linkmock,2);
link09=[ref09.galID,ref09.groupID]; %galID: CATA_Index for galaxies in the galaxy file
link09=sortrows(link09,2);
link12=[ref12.galID,ref12.groupID];
link12=sortrows(link12,2);
link15=[ref15.galID,ref15.groupID];
link15=sortrows(link15,2);

cd /home/kam/Projects/Lensing/output
%% mock halo mass
mockmass.Complexity=catmock.Gnum;  %number of different group masses
mockmass.Mmean=catmock.Gnum;   % average mass
mockmass.Gmean=catmock.Gnum;   % geometric average mass, or average in logspace
mockmass.Mmode=catmock.Gnum;   % most frequent mass
mockmass.Purity=catmock.Gnum;  % frequency of Mmode
mockmass.MBCG=mock.masstot(catmock.BCGRef); %group mass for BCG
mockmass.MIter=mock.masstot(catmock.IterRef); %group mass for BCG
GrpOffset=0;
for gid=1:numel(catmock.Gnum)
    mass=mock.masstot(linkmock(GrpOffset+(1:catmock.Mult(gid)),1));
    GrpOffset=GrpOffset+catmock.Mult(gid);
    mass=sort(mass);
    mockmass.Complexity(gid)=sum(logical(diff(mass)))+1;
    mockmass.Mmean(gid)=mean(mass);
    mockmass.Gmean(gid)=geomean(mass);
    [mockmass.Mmode(gid),mockmass.Purity(gid)]=mode(mass);
    if catmock.Mult(gid)==2
        mockmass.Mmode(gid)=max(mass);
    end
end
mockmass.Purity=mockmass.Purity./catmock.Mult;
%%
nq=max(catmock.Mult);
Q=zeros(nq,1);eQ=Q;
for i=1:nq
    Q(i)=mean(mockmass.Purity(catmock.Mult==i));
    eQ(i)=std(mockmass.Purity(catmock.Mult==i));
end
figure;
clr=colormap(hot(7));
for i=1:5
    f=catmock.Mult==i;
    if sum(f)
    [tmp,tmp1,tmp2,h]=linhist(mockmass.Purity(f),10,'stairsnorm'); hold on;
    set(h,'color',clr(i+1,:));
    end
end
myfigure;errorbar(1:nq,Q,eQ,'.');
ylabel('Average Purity');xlabel('Multiplicity');
xlim([0,nq+1]);ylim([0.4,1.1]);
print('-depsc','MockGroupQuality.eps');
%% Complexity Histograms
x=min(mockmass.Complexity)-0.5:max(mockmass.Complexity)+0.5;
myfigure;[xm,ym]=linhist(mockmass.Complexity,x);
f2=catmock.Mult==2;
f3=catmock.Mult==3;
f4=catmock.Mult>=4;
[xm2,ym2]=linhist(mockmass.Complexity(f2),x);
[xm3,ym3]=linhist(mockmass.Complexity(f3),x);
[xm4,ym4]=linhist(mockmass.Complexity(f4),x);
bar(xm,[ym;ym2;ym3;ym4]');
% set(gca,'yscale','log');
set(gca,'xtick',0:15)
xlabel('Mock Group Complexity');ylabel('Counts');
legend('All','Mult=2','Mult=3','Mult>3');
print('-depsc','MockMassComplexity.eps');
%%
% figure;loglog(mockmass.MBCG,mockmass.Mmode,'.');
% figure;loglog(mockmass.MBCG,mockmass.Mmean,'.');
f=catmock.Mult==2;

%%
f=catmock.Mult<=3;
myfigure;loglog(mockmass.Mmode(f),mockmass.Mmean(f),'g.','markersize',5);
f=catmock.Mult>3;
hold on;loglog(mockmass.Mmode(f),mockmass.Mmean(f),'ro');
legend('Mult<=3','Mult>3');legend('location','southeast');
xlabel('Mass Mode');ylabel('Mass Mean');
print('-depsc','MockMassCmpMean_Mode.eps');
%% Compare BCG mass and Mmode
fBCG=mockmass.MBCG==mockmass.Mmode;
x=1:max(catmock.Mult);
y=zeros(size(x));
for i=x
    y(i)=sum(fBCG(catmock.Mult==i))/numel(fBCG(catmock.Mult==i));
end
myfigure;plot(x,y,'o');
xlabel('Multiplicity');ylabel('BCG-Meaningful Fraction');
hold on;
plot(x,repmat(sum(fBCG)/numel(fBCG),size(x)),'-');
plot(x,Q,'r.');
l=legend('BCG','Total','Purity');
set(l,'location','southeast');
print('-depsc','MockGroupBCGRepFraction.eps');
% figure;loglog(mockmass.MBCG(~fBCG),mockmass.Mmode(~fBCG),'.')
%% Purity histogram
f=catmock.Mult>=2;
myfigure;linhist(mockmass.Purity(f),20,'stairsnorm');xlabel('Mock Group Mass Purity');ylabel('Fraction');
print('-depsc','MockMassPurity.eps');
%% Mass distributions
myfigure;loghist(mockmass.Mmean,30,'stairs','r');
hold on;loghist(mockmass.Mmode,30,'stairs','g');
loghist(mockmass.MBCG,30,'stairs','b');
legend('Mean','Mode','BCG');xlabel('Mass');ylabel('Counts');
print('-depsc','MockMassDistrMean_Mode.eps');
%% Min and Max mass for 2-member
f=catmock.Mult<=2;
myfigure;loglog(mockmass.Mmode(f),2*mockmass.Mmean(f)-mockmass.Mmode(f),'b.','markersize',3);
xlabel('Mmax');ylabel('Mmin');
title('Mult=2');
print('-depsc','MockMassMinMax.eps')
%% mass observables
f=catmock.Mult>2;
% figure;loglog(mockmass.Mmode(f),catmock.Mass(f),'.');
% figure;loglog(mockmass.Mmean(f),catmock.Mass(f),'.');
myfigure;loglog(mockmass.MBCG(f),catmock.LumMass(f),'r.','markersize',5);hold on;loglog(mockmass.MBCG(f),catmock.Mass(f),'g.','markersize',5);
loglog(mockmass.MBCG(f),sqrt(catmock.Mass(f).*catmock.LumMass(f)),'b.','markersize',5);
xlim([1e10,1e17]);ylim([1e10,1e17]);legend('Luminosity Mass','Dynamical Mass','Combined Mass');legend('location','northwest');
xlabel('$M_{BCGhost}/(M_{\odot}/h)$');ylabel('$Mass(M_{\odot}/h)$');title('$Mult>2$');
plot([1e10,1e17],[1e10,1e17],'r-');
% print('-depsc','MassObservable.eps');
f=catmock.Mult==2;
% figure;loglog(mockmass.Mmode(f),catmock.Mass(f),'.');
% figure;loglog(mockmass.Mmean(f),catmock.Mass(f),'.');
myfigure;loglog(mockmass.MBCG(f),catmock.LumMass(f),'r.','markersize',4);hold on;loglog(mockmass.MBCG(f),catmock.Mass(f),'g.','markersize',4);
% loglog(mockmass.MBCG(f),sqrt(catmock.Mass(f).*catmock.LumMass(f)),'b.','markersize',4);
xlim([1e10,1e17]);ylim([1e10,1e17]);legend('Luminosity Mass','Dynamical Mass','Combined Mass');legend('location','northwest');
xlabel('$M_{BCGhost}/(M_{\odot}/h)$');ylabel('$Mass(M_{\odot}/h)$');title('$Mult=2$');
plot([1e10,1e17],[1e10,1e17],'r-');
print('-depsc','MassObservable2.eps');
%% mass observable--best construction
mockmass.Mbest=catmock.Gnum;
mockmass.Mbest(catmock.Mult<3)=catmock.LumMass(catmock.Mult<3);
mockmass.Mbest(catmock.Mult>=3)=sqrt(catmock.LumMass(catmock.Mult>=3).*catmock.Mass(catmock.Mult>=3));
myfigure;loglog(mockmass.MBCG(f),mockmass.Mbest(f),'b.','markersize',4);hold on;
xlim([1e10,1e17]);ylim([1e10,1e17]);
xlabel('$M_{BCGhost}/(M_{\odot}/h)$');ylabel('$Mbest(M_{\odot}/h)$');
title({'$M_{best}=M_{Lum}~(Mult=2),=\sqrt(M_{Lum}M_{dyn})~(Mult>=3)$'})
plot([1e10,1e17],[1e10,1e17],'r-');
print('-depsc','MassObservable3.eps');
%%
[xx,yy,n,s]=densitygrid(mockmass.Complexity,mockmass.Purity,[16,20]);
figure;pcolor(xx,yy,log10(n+1));xlabel('Complexity');ylabel('Purity');colormap('gray');
[xx,yy,n,s]=densitygrid(catmock.Mult,mockmass.Purity,[153,20]);
figure;pcolor(xx,yy,log10(n+1));xlabel('Multiplicity');ylabel('Purity');colormap('gray');
[xx,yy,n]=gridcounts(catmock.Mult,mockmass.Complexity);
figure;pcolor(xx-0.5,yy-0.5,log10(n+1));xlabel('Multiplicity');ylabel('Complexity');colormap('gray');
print('-depsc','MockGroupComplexity-MultJoint.eps');
%%
f=catmock.Mult>1;
% figure;loglog(mockmass.Mmode(f),catmock.Mass(f),'.');
% figure;loglog(mockmass.Mmean(f),catmock.Mass(f),'.');
myfigure;loglog(mockmass.Mmode(f),catmock.LumMass(f),'r.','markersize',5);hold on;loglog(mockmass.Mmode(f),catmock.Mass(f),'g.','markersize',5);
loglog(mockmass.Mmode(f),sqrt(catmock.Mass(f).*catmock.LumMass(f)),'b.','markersize',5);
xlim([1e10,1e17]);ylim([1e10,1e17]);legend('Luminosity Mass','Dynamical Mass','Combined Mass');legend('location','northwest');
xlabel('$M_{halo}/(M_{\odot}/h)$');ylabel('$Mass(M_{\odot}/h)$');title('$Mult>2$');
plot([1e10,1e17],[1e10,1e17],'r-');
% print('-depsc','MassObservable.eps');
myfigure;loglog(mockmass.Mmode(f),mockmass.Mbest(f),'b.','markersize',5);hold on;
loglog(mockmass.Mmode(f),catmock.LumMass(f),'k.','markersize',5);
xlim([1e10,1e17]);ylim([1e10,1e17]);
xlabel('$M_{halo}/(M_{\odot}/h)$');ylabel('$Mass(M_{\odot}/h)$');title('$Mult>2$');
plot([1e10,1e17],[1e10,1e17],'r-');
% figure;loglog(catmock.LumMass(f),catmock.Mass(f),'.');
% figure;loglog(mockmass.Mmode(f),(catmock.Mass(f)+catmock.LumMass(f))/2,'b.','markersize',5);xlim([1e10,1e17]);ylim([1e10,1e17]);
% figure;loglog(mockmass.Mmean(f),catmock.TotFluxInt(f),'.');
% figure;loglog(mockmass.Mmean(f),catmock.TotFlux(f),'.');
% figure;loglog(mockmass.Mmean(f),catmock.BCGmag(f),'.');
% figure;loglog(mockmass.Mmean(f),catmock.VelDisp(f),'.');
% figure;loglog(mockmass.Mmean(f),catmock.Rad1Sig(f),'.');
%% group redshift
mockzz.Mean=catmock.Gnum;   % average mass
mockzz.Med=catmock.Gnum;   % most frequent mass
mockzz.Disp=catmock.Gnum;  % frequency of Mmode
mockzz.BCG=catmock.Gnum;
mockzz.Ext=catmock.Gnum;
GrpOffset=0;
for gid=1:numel(catmock.Gnum)
    zz=mock.zz(linkmock(GrpOffset+(1:catmock.Mult(gid)),1));
    mockzz.BCG(gid)=zz(linkmock(GrpOffset+(1:catmock.Mult(gid)),1)==catmock.BCGRef(gid));
    mockzz.Mean(gid)=mean(zz);
    mockzz.Med(gid)=median(zz);
    mockzz.Disp(gid)=std(zz);
    mockzz.Ext(gid)=max(zz)-min(zz);
    GrpOffset=GrpOffset+catmock.Mult(gid);
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
f=mockmass.Mmode>0;
f=f&catmock.Mult>=2;
myfigure;
[xx,yy,n,s]=densitygrid(log10(mockmass.Mmode(f)),catmock.MedianZ(f),[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(catmock.VelDisp'),minmax(catmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$log(M_{dyn}/(M_{\odot}/h))$');
ylabel('$z$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
% print('-depsc','mockGroupMassRedshift_3.eps');
%%
nbin=5;
%% Multiplicity fun
myfigure;
% [xm,ym,dym]=loghist(catmock.Mult,nbin,'stairs','k--');
% hold on;
% [x1,y1,dy1]=loghist(cat09.Mult,nbin,'stairs','r');
% [x2,y2,dy2]=loghist(cat12.Mult,nbin,'stairs','g');
% [x3,y3,dy3]=loghist(cat15.Mult,nbin,'stairs','b');
[x,y,dy]=loghist(cat.Mult(cat.Mult>2),nbin,'stairs','k');
xlabel('Group Multiplicity');
ylabel('dN');
hold on;
x=[3,4,5,6,7,8,9,12,18,26,41,71,220];
y=[58788,27083,14925,8744,5630,3858,6196,4427,1711,787,272,47,0];
hold on;stairs(x,y,'m');
% legend('mock','G09','G12','G15','GAMA','Jhonston','location','best');
legend('GAMA','Jhonston','location','best');
print('-depsc','GroupMultfunCmp.eps');
%%
dym=dym/sum(ym);dy1=dy1./sum(y1);dy2=dy2./sum(y2);dy3=dy3/sum(y3);dy=dy/sum(y);
myfigure;
loglog(xm,dym,'k.--');
hold on;
loglog(x1,dy1,'r.-');
loglog(x2,dy2,'g.-');
loglog(x3,dy3,'b.-');
loglog(x,dy,'k.-');
xlabel('Group Multiplicity');
ylabel('$dN/dMult/N_{tot}$');
legend('mock','G09','G12','G15','GAMA','location','best');
% print('-depsc','GroupMultfunSpec.eps');
%% Massfun
myfigure;
[xm,ym,dym]=loghist(catmock.Mass,nbin,'stairs','k--');
hold on;
[x1,y1,dy1]=loghist(cat09.Mass,nbin,'stairs','r');
[x2,y2,dy2]=loghist(cat12.Mass,nbin,'stairs','g');
[x3,y3,dy3]=loghist(cat15.Mass,nbin,'stairs','b');
[x,y,dy]=loghist([cat09.Mass;cat12.Mass;cat15.Mass],nbin,'stairs','k');
xlabel('$M_{dyn}/(M_{\odot}/h)$');
ylabel('dN');
legend('mock','G09','G12','G15','GAMA','location','best');
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
legend('mock','G09','G12','G15','GAMA','location','best');

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
[xm,ym,dym]=loghist(catmock.LumMass,nbin,'stairs','k--');
hold on;
[x1,y1,dy1]=loghist(cat09.LumMass,nbin,'stairs','r');
[x2,y2,dy2]=loghist(cat12.LumMass,nbin,'stairs','g');
[x3,y3,dy3]=loghist(cat15.LumMass,nbin,'stairs','b');
[x,y,dy]=loghist(cat.LumMass,nbin,'stairs','k');
xlabel('$M_{lum}/(M_{\odot}/h)$');
ylabel('dN');
legend('mock','G09','G12','G15','GAMA','location','best');
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
legend('mock','G09','G12','G15','GAMA','location','best');

print('-depsc','GroupLumMassfunSpec.eps');
%% Lumosity Fun
myfigure;
[xm,ym,dym]=linhist(catmock.TotFluxInt,nbin,'stairs','k--');
hold on;
[x1,y1,dy1]=linhist(cat09.TotFluxInt,nbin,'stairs','r');
[x2,y2,dy2]=linhist(cat12.TotFluxInt,nbin,'stairs','g');
[x3,y3,dy3]=linhist(cat15.TotFluxInt,nbin,'stairs','b');
[x,y,dy]=linhist(cat.TotFluxInt,nbin,'stairs','k');
set(gca,'yscale','log');
xlabel('$M_{r}$');
ylabel('dN');
legend('mock','G09','G12','G15','GAMA','location','south');
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
legend('mock','G09','G12','G15','GAMA','location','south');
print('-depsc','GroupLumfunSpec.eps');
%% Redshift Distr
myfigure;
linhist(catmock.MedianZ,nbin,'stairs','k--');
hold on;
linhist(cat09.MedianZ,nbin,'stairs','r');
linhist(cat12.MedianZ,nbin,'stairs','g');
linhist(cat15.MedianZ,nbin,'stairs','b');
linhist([cat09.MedianZ;cat12.MedianZ;cat15.MedianZ],nbin,'stairs','k');
set(gca,'yscale','log');
xlabel('$z$');
ylabel('dN');
legend('mock','G09','G12','G15','GAMA','location','best');

print('-depsc','GroupRedshift.eps');
%% Redshift Completeness
xbin=(0.05:0.1:1.05)';
myfigure;
linhist(catmock.ZComp,xbin,'stairsnorm','k--');
hold on;
linhist(cat09.ZComp,xbin,'stairsnorm','r');
linhist(cat12.ZComp,xbin,'stairsnorm','g');
linhist(cat15.ZComp,xbin,'stairsnorm','b');
linhist(cat.ZComp,xbin,'stairsnorm','k');
xlabel('Completeness');
ylabel('Fraction');
legend('mock','G09','G12','G15','GAMA','location','best');

print('-depsc','GroupZComplete.eps');
%% VD Distr
myfigure;
[xm,ym,dym]=loghist(catmock.VelDisp,nbin,'stairs','k--');
hold on;
[x1,y1,dy1]=loghist(cat09.VelDisp,nbin,'stairs','r');
[x2,y2,dy2]=loghist(cat12.VelDisp,nbin,'stairs','g');
[x3,y3,dy3]=loghist(cat15.VelDisp,nbin,'stairs','b');
[x,y,dy]=loghist(cat.VelDisp,nbin,'stairs','k');
xlabel('$\sigma/(km/s)$');
ylabel('dN');
legend('mock','G09','G12','G15','GAMA','location','best');
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
legend('mock','G09','G12','G15','GAMA','location','south');
print('-depsc','GroupVelDispSpec.eps');
%% Size Distr
myfigure;
[xm,ym,dym]=loghist(catmock.Rad50,nbin,'stairs','k--');
hold on;
[x1,y1,dy1]=loghist(cat09.Rad50,nbin,'stairs','r');
[x2,y2,dy2]=loghist(cat12.Rad50,nbin,'stairs','g');
[x3,y3,dy3]=loghist(cat15.Rad50,nbin,'stairs','b');
[x,y,dy]=loghist(cat.Rad50,nbin,'stairs','k');
xlabel('$R_{0.5}/(Mpc/h)$');
ylabel('dN');
legend('mock','G09','G12','G15','GAMA','location','best');
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
legend('mock','G09','G12','G15','GAMA','location','southwest');
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
[xx,yy,n,s]=densitygrid(catmock.VelDisp,catmock.Rad50,[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(catmock.VelDisp'),minmax(catmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$\sigma/(km/s)$');
ylabel('$R_{0.5}/(Mpc/h)$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','mockGroupSizeVD.eps');

myfigure;
[xx,yy,n,s]=densitygrid(cat.VelDisp,cat.Rad50,[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(catmock.VelDisp'),minmax(catmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$\sigma/(km/s)$');
ylabel('$R_{0.5}/(Mpc/h)$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','GAMAGroupSizeVD.eps');
%% Size+Redshift joint distr
nbin=30;
myfigure;
[xx,yy,n,s]=densitygrid(catmock.Rad50,catmock.MedianZ,[nbin,nbin]);
% pcolor(xx,yy,log10(n+1));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(catmock.VelDisp'),minmax(catmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
ylabel('$z$');
xlabel('$R_{0.5}/(Mpc/h)$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','mockGroupSizeRedshift.eps');

myfigure;
[xx,yy,n,s]=densitygrid(cat.Rad50,cat.MedianZ,[nbin,nbin]);
% pcolor(xx,yy,log10(n+1));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(catmock.VelDisp'),minmax(catmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
ylabel('$z$');
xlabel('$R_{0.5}/(Mpc/h)$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','GAMAGroupSizeRedshift.eps');
%% VD+Redshift joint distribution
nbin=30;
myfigure;
[xx,yy,n,s]=densitygrid(catmock.VelDisp,catmock.MedianZ,[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(catmock.VelDisp'),minmax(catmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$\sigma/(km/s)$');
ylabel('$z$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','mockGroupVDRedshift.eps');

myfigure;
[xx,yy,n,s]=densitygrid(cat.VelDisp,cat.MedianZ,[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(catmock.VelDisp'),minmax(catmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$\sigma/(km/s)$');
ylabel('$z$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','GAMAGroupVDRedshift.eps');
%% Mass+Redshift joint distribution
nbin=30;
f=catmock.Mass>0;
f=f&catmock.Mult>=2;
myfigure;
[xx,yy,n,s]=densitygrid(log10(catmock.Mass(f)),catmock.MedianZ(f),[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(catmock.VelDisp'),minmax(catmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$log(M_{dyn}/(M_{\odot}/h))$');
ylabel('$z$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
% print('-depsc','mockGroupMassRedshift_3.eps');

f=cat.Mass>0;
f=f&cat.Mult>=2;
myfigure;
[xx,yy,n,s]=densitygrid(log10(cat.Mass(f)),cat.MedianZ(f),[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(catmock.VelDisp'),minmax(catmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$log(M_{dyn}/(M_{\odot}/h))$');
ylabel('$z$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
% print('-depsc','GAMAGroupMassRedshift_3.eps');
%% LumMass+Redshift joint distribution
nbin=30;
f=catmock.LumMass>0;
f=f&catmock.Mult>=3;
myfigure;
[xx,yy,n,s]=densitygrid(log10(catmock.LumMass(f)),catmock.MedianZ(f),[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(catmock.VelDisp'),minmax(catmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$log(M_{Lum}/(M_{\odot}/h))$');
ylabel('$z$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','mockGroupLumMassRedshift_3.eps');

f=cat.LumMass>0;
f=f&cat.Mult>=3;
myfigure;
[xx,yy,n,s]=densitygrid(log10(cat.LumMass(f)),cat.MedianZ(f),[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(catmock.VelDisp'),minmax(catmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$log(M_{Lum}/(M_{\odot}/h))$');
ylabel('$z$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','GAMAGroupLumMassRedshift_3.eps');
%% Color-Mag joint distribution
nbin=30;
% myfigure;
% [xx,yy,n,s]=densitygrid(mock.AB_u-mock.AB_r,mock.AB_r,[nbin,nbin]);
% % pcolor(xx,yy,log10(1+n));hold on;
% contourf(xx,yy,log10((n+1)))
% % imagesc(minmax(catmock.VelDisp'),minmax(catmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
% colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
% xlabel('$U_{AB}-R_{AB}$');
% ylabel('$R_{AB}$');
% title([num2str(nbin),'$\times$',num2str(nbin)]);
% print('-depsc','mockColorMag.eps');

myfigure;
[xx,yy,n,s]=densitygrid(gal.AB_u-gal.AB_r,gal.AB_r,[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(catmock.VelDisp'),minmax(catmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$U_{AB}-R_{AB}$');
ylabel('$R_{AB}$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','GAMAColorMag.eps');
%% Mass-Mult
f=catmock.Mass>0;
myfigure;
loglog(catmock.Mass(f),catmock.Mult(f),'.','markersize',5);
xlabel('$M_{dyn}/(M_{\odot}/h)$');
ylabel('$Multiplicity$');
print('-depsc','mockGroupMassMult.eps');

f=cat.Mass>0;
myfigure;
loglog(cat.Mass(f),cat.Mult(f),'.','markersize',5);
xlabel('$M_{dyn}/(M_{\odot}/h)$');
ylabel('$Multiplicity$');
print('-depsc','GAMAGroupMassMult.eps');
%%
myfigure;
loglog(catmock.VelDisp,catmock.Mass,'.','markersize',3);
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
