%%
% view GAMA catalogue G3Cv4
%%
load_G3Cv4;
[grp09,grp12,grp15]=split_mockcat(grp);
%%
f=gal.LumMass>0&~isnan(gal.logSFR);
f1=gal.IsIterCen==0&f;
f2=gal.IsIterCen>0&f;
myfigure;loglog(gal.LumMass(f1),10.^gal.logSFR(f1),'r.');
hold on;
loglog(gal.LumMass(f2),10.^gal.logSFR(f2),'g.');
ylim([1e-3,1e3]);
xlabel('LumMass/Msun');ylabel('SFR/(Msun/yr)');
legend('Central','Sat');
print('-depsc','SFR_Mhost.eps');
%%
f=gal.LumMass>0&~isnan(gal.logSSFR);
f1=gal.IsIterCen==0&f;
f2=gal.IsIterCen>0&f;
myfigure;loglog(gal.LumMass(f1),10.^gal.logSSFR(f1),'r.');
hold on;
loglog(gal.LumMass(f2),10.^gal.logSSFR(f2),'g.');
ylim([1e-15,1e-5]);
xlabel('LumMass/Msun');ylabel('SFR/(Msun/yr)');
legend('Central','Sat');
print('-depsc','SSFR_Mhost.eps');
%%
figure;plot([C1(:),C2(:),C3(:)]);
f=halomock.Mult>1;
myfigure;loglog(halomock.HaloMass(f),halomock.MassProxy(f),'.');
hold on;
loglog(halomock.HaloMass(f),halomock.LumMass(f),'r.');
plot([1e10,1e15],[1e10,1e15]);
% loglog(halomock.HaloMass,halomock.HybMass,'g.');
xlabel('Halo Mass [M$_\odot/h$]');
ylabel('Mass Proxy [M$_\odot/h$]');
l=legend('Dynamical','Luminosity');set(l,'location','southeast');
print('-depsc','MassHaloMock.eps');
%%
f=grpmock.Mult>1;
myfigure;loglog(mockmass.MIter(f),grpmock.MassProxy(f),'.');
hold on;
loglog(mockmass.MIter(f),grpmock.LumMass(f),'r.');
plot([1e10,1e15],[1e10,1e15]);
% loglog(halomock.HaloMass,halomock.HybMass,'g.');
xlabel('Group Mass [M$_\odot/h$]');
ylabel('Mass Proxy [M$_\odot/h$]');
l=legend('Dynamical','Luminosity');set(l,'location','southeast');
print('-depsc','MassGroupMock.eps');
%%
f=grp.Mult>1;
myfigure;loglog(grp.MassProxy(f),grp.LumMass(f),'.');
hold on;
plot([1e9,1e17],[1e9,1e17]);
xlabel('Dynamical Mass [M$_\odot/h$]');
ylabel('Luminosity Mass [M$_\odot/h$]');
xlim([1e9,1e17]);
ylim([1e9,1e17]);
print('-depsc','MassObs.eps');
%%
f=grpmock.Mult>1;
figure;
loglog(grpmock.VelDisp(f),grpmock.MassProxy(f),'.');
xlabel('\sigma_v');ylabel('DynMass');
print('-depsc','DynMass_VelDisp.eps');
figure;
loglog(grpmock.Rad50(f),grpmock.MassProxy(f),'.');
xlabel('Rad50 (Mpc/h)');ylabel('DynMass');
xlim([1e-11,1e2]);
print('-depsc','DynMass_Rad50.eps');
%%
figure;loglog(grpmock.LumMass,grpmock.MassProxy,'.');
hold on;
plot([1e10,1e15],[1e10,1e15]);

figure;loglog(grp.LumMass,grp.MassProxy,'.');
hold on;
plot([1e10,1e15],[1e10,1e15]);

figure;loglog(grp.LumMass,grp.TotFluxProxy,'.');
figure;semilogx(grp.TotFluxProxy,grp.LumMass./grp.TotFluxProxy,'.')
%% Choose one from below
dir='/home/kam/Projects/Lensing/data/2011_02';
Tz0=importdata([dir,'/Tz0.new'],' ',1);
galmock.HaloMass=Tz0.data(:,8);
galmock.Mag=Tz0.data(:,5);
%%  Choose this one or above.
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
    mag=galmock.Mag(linkmock(GrpOffset+(1:grpmock.Mult(gid)),1));
    ids=galmock.HaloID(linkmock(GrpOffset+(1:grpmock.Mult(gid)),1)); %hostids for member galaxies in the group
    ids=[ids,galmock.GalID(linkmock(GrpOffset+(1:grpmock.Mult(gid)),1))];%[HaloID, GalID]
    GrpOffset=GrpOffset+grpmock.Mult(gid);
    ids=sortrows(ids);
    mockmass.Complexity(gid)=sum(logical(diff(ids(:,1))))+1;
    mockmass.Mmean(gid)=mean(mass);
    mockmass.Gmean(gid)=geomean(mass(logical(mass)));  %geomean of non-zero elements
    [mockmass.Mmode(gid),mockmass.Purity(gid)]=mode(ids(:,1));
    mid=ids(ids(:,1)==mockmass.Mmode(gid),2);mid=mid(1); %take GalID of the mode galaxy
    mockmass.Mmode(gid)=galmock.HaloMass(mid);
    if mockmass.Purity(gid)==1
%         mockmass.Mmode(gid)=max(mass);
        [~,minmagid]=min(mag);
        mockmass.Mmode(gid)=mass(minmagid);
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
figure;loglog(mockmass.MIter,mockmass.Mmode,'.');
figure;loglog(mockmass.MIter,mockmass.Mmean,'.');
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
print('-depsc','MockGroupBCGRepFraction.eps');
figure;loglog(mockmass.MBCG(~fBCG),mockmass.Mmode(~fBCG),'.')
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
print('-depsc','MockGroupIterRepFraction.eps');
% figure;loglog(mockmass.MIter(~fIter),mockmass.Mmode(~fIter),'.')
%% Purity histogram
f=grpmock.Mult>=2;
myfigure;linhist(mockmass.Purity(f),20,'stairsnorm');xlabel('galmock Group Mass Purity');ylabel('Fraction');
print('-depsc','MockMassPurity.eps');
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
loglog(mockmass.MIter(f),grpmock.MassProxy(f),'g.','markersize',5);hold on;
loglog(mockmass.MIter(f),grpmock.LumMass(f),'r.','markersize',5);
% loglog(mockmass.MBCG(f),sqrt(grpmock.MassProxy(f).*grpmock.LumMass(f)),'b.','markersize',5);
% xlim([1e10,1e17]);ylim([1e10,1e17]);
legend('Dynamical Mass','Luminosity Mass','Hybrid Mass');legend('location','northwest');
xlabel('$M_{Iterhost}/(M_{\odot}/h)$');ylabel('$Mass(M_{\odot}/h)$');title('$Mult>2$');
plot([1e10,1e17],[1e10,1e17],'r-');
% print('-depsc','MassObservable.eps');
%%
nbin=30;
myfigure;
[xx,yy,n,s]=densitygrid(log10(mockmass.MIter),log10(grpmock.LumMass),[nbin,nbin],[10,17],[10,17]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(grpmock.VelDisp'),minmax(grpmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$IterMass$');
ylabel('$ObsMass$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
nbin=30;
myfigure;
[xx,yy,n,s]=densitygrid(log10(mockmass.MIter),log10(grpmock.MassProxy),[nbin,nbin],[10,17],[10,17]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(grpmock.VelDisp'),minmax(grpmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$IterMass$');
ylabel('$ObsMass$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
nbin=30;
myfigure;
[xx,yy,n,s]=densitygrid(log10(mockmass.MIter),log10(sqrt(grpmock.MassProxy.*grpmock.LumMass)),[nbin,nbin],[10,17],[10,17]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(grpmock.VelDisp'),minmax(grpmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$IterMass$');
ylabel('$ObsMass$');
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
print('-depsc','MassObservable2.eps');
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
myfigure;loglog(mockmass.Mmode(f),grpmock.LumMass(f),'r.','markersize',5);hold on;loglog(mockmass.Mmode(f),grpmock.MassProxy(f),'g.','markersize',5);
loglog(mockmass.Mmode(f),sqrt(grpmock.MassProxy(f).*grpmock.LumMass(f)),'b.','markersize',5);
xlim([1e10,1e17]);ylim([1e10,1e17]);legend('Luminosity Mass','Dynamical Mass','Combined Mass');legend('location','northwest');
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
% print('-depsc','mockGroupMassRedshift.eps');
%%
nbin=5;
%% Multiplicity fun
myfigure;
% [xm,ym,dym]=loghist(grpmock.Mult,nbin,'stairs','k--');
% hold on;
% [x1,y1,dy1]=loghist(grp09.Mult,nbin,'stairs','r');
% [x2,y2,dy2]=loghist(grp12.Mult,nbin,'stairs','g');
% [x3,y3,dy3]=loghist(grp15.Mult,nbin,'stairs','b');
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
%% Massfun show
nbin=logspace(9,17,17);
f=grp.Mult>2;
myfigure;
[xm,ym,dym]=loghist(grp.MassProxy(f),nbin);
hold on;
[x,y,dy]=loghist(grp.LumMass(f),nbin);
ym(ym==0)=0.1;
y(y==0)=0.1;
stairs(nbin,[ym,0.1],'k-');
stairs(nbin,[y,0.1],'k--');
box on;
xlim([1e9,1e17]);ylim([1,1e4]);
set(gca,'xtick',10.^(9:17));
set(gca,'xscale','log');
set(gca,'yscale','log');
xlabel('Mass[$M_\odot/h$]');
ylabel('Counts');
legend('Dynamical Mass', 'Luminosity Mass');
print('-depsc','GroupMassFun.eps');
%% Massfun
myfigure;
ax=GAMAloghist(grpmock,grp,20,'LumMass');
% set(ax(1),'ylim',[10,2e4]);
% set(ax(2),'ylim',[1e-22,1e-4]);
set(get(ax(2),'xlabel'),'string','LumMass/Msun');
set(get(ax(2),'ylabel'),'string','dN/dM/Ntot');
print('-depsc','GroupLumMassfun.eps');

myfigure;
ax=GAMAloghist(grpmock,grp,20,'MassProxy');
set(get(ax(2),'xlabel'),'string','$DynMass/Msun$');
set(get(ax(2),'ylabel'),'string','dN/dM/Ntot');
print('-depsc','GroupDynMassfun.eps');

myfigure;
ax=GAMAloghist(grpmock,grp,20,'TotFluxProxy');
set(get(ax(2),'xlabel'),'string','Luminosity/Lsun');
set(get(ax(2),'ylabel'),'string','dN/dL/Ntot');
print('-depsc','GroupLumfun.eps');

myfigure;
ax=GAMAloghist(grpmock,grp,20,'VelDisp');
set(get(ax(2),'xlabel'),'string','$\sigma_v/(km/s)$');
set(get(ax(2),'ylabel'),'string','$dN/d\sigma/Ntot$');
print('-depsc','GroupVelDispFun.eps');

%% Redshift Distr
myfigure;
linhist(grpmock.Zfof,nbin,'stairs','k--');
hold on;
linhist(grp09.Zfof,nbin,'stairs','r');
linhist(grp12.Zfof,nbin,'stairs','g');
linhist(grp15.Zfof,nbin,'stairs','b');
linhist([grp09.Zfof;grp12.Zfof;grp15.Zfof],nbin,'stairs','k');
set(gca,'yscale','log');
xlabel('$z$');
ylabel('dN');
legend('galmock','G09','G12','G15','GAMA','location','best');

print('-depsc','GroupRedshift.eps');
%% Redshift Completeness
xbin=(0.05:0.1:1.05)';
myfigure;
linhist(grpmock.Zcomp,xbin,'stairsnorm','k--');
hold on;
linhist(grp09.Zcomp,xbin,'stairsnorm','r');
linhist(grp12.Zcomp,xbin,'stairsnorm','g');
linhist(grp15.Zcomp,xbin,'stairsnorm','b');
linhist(grp.Zcomp,xbin,'stairsnorm','k');
xlabel('Completeness');
ylabel('Fraction');
legend('Mock','G09','G12','G15','GAMA','location','best');

print('-depsc','GroupZComplete.eps');
%% Size Distr
myfigure;
[xm,ym,dym]=loghist(grpmock.Rad50,nbin,'stairs','k--');
hold on;
[x1,y1,dy1]=loghist(grp09.Rad50,nbin,'stairs','r');
[x2,y2,dy2]=loghist(grp12.Rad50,nbin,'stairs','g');
[x3,y3,dy3]=loghist(grp15.Rad50,nbin,'stairs','b');
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
%%
f=grp.MassProxy>0;
f=f&grp.Mult>=2;
myfigure;
[xx,yy,n,s]=densitygrid(log10(grp.MassProxy(f)),grp.Zfof(f),[nbin,nbin]);
% pcolor(xx,yy,log10(1+n));hold on;
contourf(xx,yy,log10((n+1)))
% imagesc(minmax(grpmock.VelDisp'),minmax(grpmock.Rad50'),log10(n+0.1));set(gca,'ydir','normal');
colorbar;l=legend('$log(1+dN)$');set(l,'interpreter','latex');
xlabel('$log(M_{dyn}/(M_{\odot}/h))$');
ylabel('$z$');
title([num2str(nbin),'$\times$',num2str(nbin)]);
print('-depsc','GAMAGroupMassRedshift.eps');
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
%%
f=grp.LumMass>0;
f=f&grp.Mult>=3;
myfigure;
[xx,yy,n,s]=densitygrid(log10(grp.LumMass(f)),grp.Zfof(f),[nbin,nbin]);
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
%%
f=grp.LumMass>0;
myfigure;
loglog(grp.LumMass(f),grp.Mult(f),'.','markersize',5);
xlabel('$M_{lum}/(M_{\odot}/h)$');
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
plot(grp.CenRA,grp.CenDEC,'.')
