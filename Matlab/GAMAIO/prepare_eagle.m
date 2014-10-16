file='/work/Projects/Lensing/data/Eagle.hdf5';
Mhalo=h5read(file,'/FOF/GroupMass');
CenSub=h5read(file,'/FOF/FirstSubhaloID');
CenSub=CenSub+1;
Msub=h5read(file,'/Subhalo/Mass');
Msub=Msub(CenSub);
Msubtype=h5read(file,'/Subhalo/MassType')';
Msubtype=Msubtype(CenSub,:);
r=[1,3,5,10,20,30,40,50,70,100];
for i=1:numel(r)
    Map(:,:,i)=h5read(file,['/Subhalo/ApertureMeasurements/Mass/',num2str(r(i),'%03d'),'kpc'])';
    SFRs(:,i)=h5read(file,['/Subhalo/ApertureMeasurements/SFR/',num2str(r(i),'%03d'),'kpc']);
end    
SFRs=SFRs(CenSub,:);
%%
Mstars=squeeze(Map(CenSub,5,:));
Mgases=squeeze(Map(CenSub,1,:));
Mdm=squeeze(Map(CenSub,2,:));
%%
Mstarhalo=Mdm.*repmat(Msubtype(:,5)./Msubtype(:,2),1,numel(r));
[~,icut]=min(abs(log(diff(Mstars'))-log(diff(Mstarhalo'))));
ipick=sub2ind(size(Mstars),1:numel(icut),icut+1)';
Mgal=Mstars(ipick);
SFR=SFRs(ipick);
Mgalmax=Mstars(:,end);
SFRmax=SFRs(:,end);
%%
SFR=SFR*1e9*0.7.^2; %Msun/h^2/Gyr;
Mgal=Mgal*1e10*0.7; %Msun/h^2;
Mhalo=Mhalo*1e10; %Msun/h
Mgalmax=Mgalmax*1e10*0.7;
SFRmax=SFRmax*1e9*0.7.^2;
Eagle.Mgal=Mgal;
Eagle.Mhalo=Mhalo;
Eagle.SFR=SFR;
Eagle.SSFR=SFR./Mgal;
Eagle.Mgalmax=Mgalmax;
Eagle.SFRmax=SFRmax;
Eagle.SSFRmax=SFRmax./Mgalmax;
save /work/Projects/Lensing/data/Eagle.mat Eagle
%%
haloid=300;
figure;
loglog(r(1:end-1),diff(Mstars(haloid,:))./diff(r.^3),'r')
hold on;
loglog(r(1:end-1),diff(Mgases(haloid,:))./diff(r.^3),'g')
loglog(r(1:end-1),diff(Mdm(haloid,:))./diff(r.^3)*Msubtype(haloid,5)/Msubtype(haloid,2),'b')
% yscale('linear')
plot([r(icut(haloid)),r(icut(haloid))],ylim(),'k--');
%%
figure;
loglog(Mgal,Mhalo,'.');
%%
figure;
f=Eagle.Mgal>0.5e10;
linhist(log10(Eagle.SSFRmax(f)),linspace(-5,1),'stairs');
hold on;
linhist(log10(Eagle.SSFR(f)),linspace(-5,1),'stairs','r');
%%
figure;
loglog(Mgal,SFR./Mgal,'r.');hold on;
loglog(Mgal,SFRmax./Mgalmax,'g.');
xlim([5e9,1e12]);
%%
myfigure;
xbin=8:0.5:12;
x=log10(Mgal);
y=log10(Mhalo);
f=SFR./Mgal>10^-2;
% plot(x(f),y(f),'b.');hold on;
[xmed,ymed,yl]=skeleton(x(f),y(f),xbin,0.68);
h2=plot(xmed,ymed,'-','linewidth',2,'color','g');hold on;
plot(xmed,yl,'g--');
f=SFR./Mgal<10^-2;
% plot(x(f),y(f),'r.');
[xmed,ymed,yl]=skeleton(x(f),y(f),xbin,0.68);
h2=plot(xmed,ymed,'-','linewidth',2,'color','m');hold on;
plot(xmed,yl,'m--');
xlim([8,12]);
xlabel('$\log(M_\star$[M$_\odot/h^2$])');
ylabel('$\log(M_h$[M$_\odot/h$])');
% printf('-depsc','/work/Projects/Lensing/output/paper/paperII/extra/Eagle.eps');