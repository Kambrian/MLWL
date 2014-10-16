% grp.SMpeak=grp.IterCenSMsps;
% grp.SMpeakRef=grp.IterCenRef;
% for i=1:numel(grp.GroupID)
%     grpid=grp.GroupID(i);
%     [grp.SMpeak(i),j]=max(gal.SMsps(gal.GroupID==grpid));
%     tmp=find(gal.GroupID==grpid);
%     grp.SMpeakRef(i)=tmp(j);
% end
% flag=gal.FlagSFR>0&gal.SSFR>10^-1.5;
% grp.FlagHighSSFR=zeros(size(grp.GroupID));
% for i=1:numel(grp.GroupID)
% grp.FlagHighSSFR(i)=flag(grp.IterCenRef(i)==gal.GalID);
% end
% grp.FlagHighSSFR=logical(grp.FlagHighSSFR);
% SMdiff=log10(grp.SMpeak./grp.IterCenSMsps);
%---done---
load G3Cv4up8/G3Cv4up8.mat
%% displacement
grp.SMpeakRA=gal.RA(grp.SMpeakRef);
grp.SMpeakDEC=gal.DEC(grp.SMpeakRef);
DIterPeak=ccdist([grp.IterCenRA,grp.IterCenDEC],[grp.SMpeakRA,grp.SMpeakDEC]);
DIterPeak=acos(DIterPeak).*AD_dist_flat(0.3,0,grp.Zfof).*(1+grp.Zfof); %comoving seperation
d=DIterPeak./grp.Rad100;
figure;
linhist(d(grp.FlagHighSSFR&grp.Mult>5&log(grp.SMpeak./grp.IterCenSMsps)>0.6),0:0.02:1.2,'stairsnorm')
%% mocks
% grpmock.SMpeak=grpmock.MstarIter;
% for i=1:numel(grpmock.Gnum)
%     grpid=grpmock.Gnum(i);
%     [grpmock.SMpeak(i),j]=max(galmock.Mstar(galmock.groupID==grpid));
%     tmp=find(galmock.groupID==grpid);
%     grpmock.SMpeakRef(i)=tmp(j);
% end
% flag=galmock.SSFR>10^-1.5;
% grpmock.FlagHighSSFR=flag(grpmock.IterRef);
% SMdiffmock=log10(grpmock.SMpeak./grpmock.MstarIter);
%--done----
load G3Cv4up8/mockcat_1.mat
%% displacement
grpmock.SMpeakRA=galmock.RA(grpmock.SMpeakRef);
grpmock.SMpeakDEC=galmock.DEC(grpmock.SMpeakRef);
DIterPeakmock=ccdist([grpmock.IterCenRA,grpmock.IterCenDEC],[grpmock.SMpeakRA,grpmock.SMpeakDEC]);
DIterPeakmock=acos(DIterPeakmock).*AD_dist_flat(0.3,0,grpmock.Zfof).*(1+grpmock.Zfof); %comoving seperation
d=DIterPeakmock./grpmock.Rad100;
linhist(d(grpmock.FlagHighSSFR&grpmock.Mult==6),0:0.02:1.2,'stairsnorm')
%%
% figure;
list=find(log10(grp.SMpeak./grp.IterCenSMsps)>0.5);
for j=15
    i=list(j);
    cla;
    r=grp.Rad100(i)/AD_dist_flat(0.3,0,grp.Zfof(i))/(1+grp.Zfof(i))*180/pi;
    plot_circle([grp.IterCenRA(i),grp.IterCenDEC(i)],r)
    hold on;
    f=gal.GroupID==grp.GroupID(i);
    plot(gal.RA(f),gal.DEC(f),'.');
    plot(grp.IterCenRA(i),grp.IterCenDEC(i),'rp');
    plot(grp.SMpeakRA(i),grp.SMpeakDEC(i),'gs');
    text(grp.IterCenRA(i)+0.1*r,grp.IterCenDEC(i),num2str(gal.Rmodel_abs(gal.GalID==grp.IterCenRef(i)),'%.1f'));
    text(grp.SMpeakRA(i)+0.1*r,grp.SMpeakDEC(i),num2str(gal.Rmodel_abs(grp.SMpeakRef(i)),'%.1f'));
end
%%
sym='xp';
list=find(log10(grpmock.SMpeak./grpmock.MstarIter)>0.5);
for j=10
    i=list(j);
    cla;
    r=grpmock.Rad100(i)/AD_dist_flat(0.3,0,grpmock.Zfof(i))/(1+grpmock.Zfof(i))*180/pi;
    plot_circle([grpmock.IterCenRA(i),grpmock.IterCenDEC(i)],r);
    hold on;
    f=galmock.groupID==grpmock.Gnum(i);
    plot(galmock.RA(f),galmock.DEC(f),'.');
    plot(grpmock.IterCenRA(i),grpmock.IterCenDEC(i),['r',sym(galmock.is_central(grpmock.IterRef(i))+1)],'markersize',15);
    plot(grpmock.SMpeakRA(i),grpmock.SMpeakDEC(i),['g',sym(galmock.is_central(grpmock.SMpeakRef(i))+1)],'markersize',10);
    text(grpmock.IterCenRA(i)+0.1*r,grpmock.IterCenDEC(i),num2str(galmock.Rmodel_abs(grpmock.IterRef(i)),'%.1f'));
    text(grpmock.SMpeakRA(i)+0.1*r,grpmock.SMpeakDEC(i),num2str(galmock.Rmodel_abs(grpmock.SMpeakRef(i)),'%.1f'));
end
%% Distribution of mass offset
MultMin=1;Mmin=0e10;
myfigure;
x=0:0.01:2;
[xm,ym]=linhist(SMdiff(grp.Mult>MultMin&grp.IterCenSMsps>Mmin),x);
y=cumsum(ym(end:-1:1));
y=y(end:-1:1)/y(end);
plot(x(1:end-1),y,'k-');hold on;
[xm,ym]=linhist(SMdiff(grp.FlagHighSSFR&grp.Mult>MultMin&grp.IterCenSMsps>Mmin),x);
y=cumsum(ym(end:-1:1));
y=y(end:-1:1)/y(end);
plot(x(1:end-1),y,'b-');hold on;
%mocks
[xm,ym]=linhist(SMdiffmock(grpmock.Mult>MultMin&grpmock.MstarIter>Mmin),x);
y=cumsum(ym(end:-1:1));
y=y(end:-1:1)/y(end);
plot(x(1:end-1),y,'k--');hold on;
[xm,ym]=linhist(SMdiffmock(grpmock.FlagHighSSFR&grpmock.Mult>MultMin&grpmock.MstarIter>Mmin),x);
y=cumsum(ym(end:-1:1));
y=y(end:-1:1)/y(end);
plot(x(1:end-1),y,'b--');hold on;
plot([0,1],[0.75,0.75],'k:');
yscale('log')
xlim([0,1]);
ylim([5e-3,1]);
% l=legend(['All Groups($N>',num2str(MultMin),',M_\star>10^{10}$)'],'Central-Active Groups','Mocks');set(l,'interpreter','latex');
l=legend(['All Groups($N>',num2str(MultMin),'$)'],'Central-Active Groups','Mocks');set(l,'interpreter','latex');
xlabel('$\log(M_{\star,max})-\log(M_{\star,cen})$');
ylabel('Fraction$(>\Delta \log M)$');
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/StellarMassOffset.eps');
%% Real central fractions
myfigure;
x=0:0.2:2;
y=x;
y2=x;
for i=1:numel(x)-1
    f=SMdiffmock>=x(i);%&SMdiffmock<x(i+1);
    y(i)=sum(galmock.is_central(grpmock.IterRef(f))==0)/sum(f); %numel(grpmock.IterRef);
    f=f&grpmock.FlagHighSSFR>0;
    y2(i)=sum(galmock.is_central(grpmock.IterRef(f))==0)/sum(f); %sum(grpmock.FlagHighSSFR>0)
end
y(end)=nan;y2(end)=nan;
plot(x,y,'k-');
hold on;
plot(x,y2,'b-');
xlabel('$\log(M_{\star,max})-\log(M_{\star,cen})$');
ylabel('$P_{sat}(\geq\Delta\log M)$');
legend('Mock Central','Mock Active Central');
yscale('log');
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/FakeCentralFraction_MassOffset.eps');
%%
myfigure;
x=1:10;
y=x;
y2=x;
for i=1:numel(x)
    f=grpmock.Mult>x(i);
    y(i)=sum(galmock.is_central(grpmock.IterRef(f))==0)./sum(f);
    f=grpmock.FlagHighSSFR>0&grpmock.Mult>x(i);
    y2(i)=sum(galmock.is_central(grpmock.IterRef(f))==0)./sum(f);
end
plot(x,y,'k-');
hold on;
plot(x,y2,'b-');
xlabel('$N_{min}$');
ylabel('$P_{sat}(N>N_{min}$)');
legend('Mock Central','Mock Active Central');
yscale('log');
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/SatContamination/FakeCentralFraction_MultMin.eps');
%% measure SM-HM centered on most-massive galaxy in the group !

%%
myfigure;
frac=[62,19,12,7;
    62,19,11,8;
    61,3,31,5;
    61,3,26,10];
bar(frac,'stack');
set(gca,'xticklabel',{'All-Iter','All-Peak','Act-Iter','Act-Peak'});
ylabel('Percents');
text(0.6,30,'SameCentral','color','r');
text(0.6,70,'SameSat','color','r');
text(0.6,85,'Cen','color','g');
text(0.6,95,'Sat','color','g');
print('-depsc','/work/Projects/Lensing/outputv4/paper/paperII/extra/SatContamination/CentralComposition.eps');