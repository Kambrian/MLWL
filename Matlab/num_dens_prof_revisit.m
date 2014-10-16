x=0.1:0.01:10;
% y=normpdf(x,0,0.1);
% plot(x,y,'.')
%%
Hz=@(z,OmegaM) 100*sqrt(OmegaM .*(1+z).^3+(1-OmegaM));
plot(x,(Hz(x,0.3)/3e5))
% hold on;
% plot(x,3000*exp(-0.56*x)./(1./(Hz(x,0.3)/3e5))-1,'r')
%%

%%
plot(x,1.16e7/sqrt(2*pi)/5.98*exp(-(x+19.9).^2/2/5.98^2),'k')
s1=sqrt(5.98^2+0.1^2)
plot(x,1.16e7/sqrt(2*pi)/s1*exp(-(x+19.9).^2/2/s1^2),'o')
%%
cftool(x,1./(Hz(x,0.3)/3e5))
%%
dz=0.1;
% zc=x-dz;
alpha=8.8;
sigh=5.98;
sigz=0.1;
zh=-19.9;
A=1.16e7;
beta=1/A./exp(sigh^2*alpha^2/2+alpha*(zc-zh))*erfc((alpha*sigz^2+dz)/sqrt(2)/sigz)./erfc((alpha*(sigz^2+sigh^2)+dz+zc-zh)/sqrt(2*(sigz^2+sigh^2)));
% plot(zc,beta)
%%
alpha=8.0;
beta=@(sig) (alpha+0.56)/2*exp(alpha^2*sig.^2/2+alpha*0.1).*erfc((alpha*sig.^2+0.1)./(sqrt(2)*sig));
sigz=@(z) 1.6./(z.^2-14.5*z+21.7);
% sigz=@(z) spline(xm,ysig,z);
%%
H=100;
OmegaM=0.3;
h=@(z) H*sqrt(OmegaM .*(1+z).^3+(1-OmegaM));
c=3e5;
G=43007.1;
f=@(x) (2./x.*atan(sqrt(1-x.^2)./x)-2*sqrt(1-x.^2));
dz=0.1;
gamma=@(zh) quad(@(z) normpdf(z,zh,sigz(zh)).*source_zdistr(z).*h(z)/c,zh+dz,1.5)./quad(@(z) source_zdistr(z),zh+dz,1.5);
gammaz=@(zh,dz) quad(@(z) normpdf(z,zh,sigz(zh)).*source_zdistr(z).*h(z)/c,zh+dz,1.5)./quad(@(z) source_zdistr(z),zh+dz,1.5);
z=0:0.1:0.5;
y=z;
for i=1:numel(z)
    y(i)=gamma(z(i));
end
    plot(z,beta(sigz(z)).*h(z+0.1)/c,'r',z,y,'g')
%%
% cd /work/Projects/Lensing/outputv4/data/
cd /mnt/charon/Lensing/output/
colors='rgbk';
%------------
Ac=-1.2;An=20.7;Az=2.3;
Bc=0.94;Bn=-0.67;Bz=0.16;       
grp=loadGAMAcsv('/work/Projects/Lensing/data/G3Cv4/group/G3CFoFGroup194v04.dat',4);
grp.MassProxy=grp.MassProxy.*(Ac+An./sqrt(grp.Mult)+Az./sqrt(grp.Zfof));
grp.TotFluxProxy=grp.TotFluxProxy.*(Bc+Bn./sqrt(grp.Mult)+Bz./sqrt(grp.Zfof));
[grp.LumMass,C3,A3]=luminosity_mass(grp);
ngrps=numel(grp.LumMass);
mmin=[0,1.191e13,3.33e13,1e14];
mmax=[1.191e13,3.33e13,inf,inf];
% alphas=[0.5,0.5,0.47,0.42];
alphas=[0.45,0.45,0.45,0.45];
%-------------
dz=0.01;
for i=1:4
file=['WL_L',num2str(i),'.DZ0.hdf5']
% file=['WL_L4.hdf5'];
m=h5read(file,'/predict1/Mmean');
z=h5read(file,'/predict1/z');
r=h5read(file,'/shear/seperation');
n=h5read(file,'/shear/numpair');
n=double(n);
nr=h5read(file,'/rand/numpair');
nerr=h5read(file,'/rand/numpair_err');
% n=h5read(file,'/shear/weight');
% nr=h5read(file,'/rand/weight');
% nerr=h5read(file,'/rand/weight_err');
nmock=h5read(file,'/rand/numrand');
nmock=double(nmock);
rv=comoving_200b_radius(m,0.3)/1000;
% errorbar(r,n./r.^2,sqrt(n)./r.^2,[colors(i),'-']);hold on;
% plot(r,nr./r.^2,[colors(i),'--']);
rat=n./(nr.*nmock/max(nmock));
% rat_err=sqrt((1./n+(nerr./nr).^2)).*rat;
rat_err=nerr./nr;
errorbar(r/rv,rat-1,rat_err*2,[colors(i),'o']);
hold on;
% plot(r/rv,1+gammaz(z,dz).*nfw_surf_overdensity(r/(1+z),m,z,2),[colors(i),'-'])
plot(r/rv,gammaz(z,dz).*sub_surf_overdensity(r/(1+z),m,z,2,1/alphas(i)),[colors(i),'--']) %r_{-2}=0.6R200c~[0.4~0.5]R200b
% Rv=(G*m/OmegaM).^(1/3)/1000;
% req=sqrt(200/3)*Rv;
% plot(r,gamma(z).*f(r/req)*req,'r-');
%----------------
% s1=zeros(size(r));
% s2=0;
% for gid=1:ngrps
%     if grp.Mult(gid)<=2, continue;end
%     if grp.LumMass(gid)<mmin(i)||grp.LumMass(gid)>mmax(i), continue;end
%     zh=grp.Zfof(gid);
%     m=grp.LumMass(gid)/1e10;
% %     m=grp.MassProxy(gid)/1e10;
%     if isnan(m)||m==0, continue;end
% %     s1=s1+quad(@(z) normpdf(z,zh,sigz(zh)).*source_zdistr(z).*h(z)/c,zh+dz,1.5).*nfw_surf_overdensity(r/(1+zh),m,zh,1);
%     s1=s1+quad(@(z) normpdf(z,zh,sigz(zh)).*source_zdistr(z).*h(z)/c,zh+dz,1.5).*sub_surf_overdensity(r/(1+zh),m,zh,2,1/alphas(i));
% %     s1=s1+beta(sigz(z))*exp(-alpha*(z+0.1)).*nfw_surf_overdensity(r/(1+z),m,z,1);
% %     Rv=(G*m/OmegaM).^(1/3)/1000;
% %     req=sqrt(200/3)*Rv;
% %     s1=s1+beta(sigz(z))*exp(-alpha*(z+0.1)).*f(r/req)*req;
%     s2=s2+quad(@(z) source_zdistr(z),zh+dz,1.5);
% end
% plot(r,1+s1./s2,[colors(i),':']);
% plot(r,(n./(nr.*nmock/max(nmock))-1)./(s1./s2),[colors(i),'o-']);
%--------------
end
set(gca,'xscale','log')
% set(gca,'yscale','log')
%%
id=4;
bid=num2str(id);
files={['WL_L',bid,'.DZ0.hdf5'],['WL_L',bid,'.hdf5'],['WL_L',bid,'.DZ2.hdf5'],['WL_L',bid,'.DZ3.hdf5']};
dzs=[0.01,0.1,0.2,0.3];
alphas=[0.45,0.45,0.45,0.45];
% file=files{4};
% n=h5read(file,'/shear/numpair');n=double(n);nmock=h5read(file,'/rand/numrand');nmock=double(nmock);nr=h5read(file,'/rand/numpair');
% rat_res=n./(nr.*nmock/max(nmock));
for i=1:numel(files)
file=files{i}
dz=dzs(i);
% file=['WL_L4.hdf5'];
m=h5read(file,'/predict1/Mmean');
z=h5read(file,'/predict1/z');
r=h5read(file,'/shear/seperation');
n=h5read(file,'/shear/numpair');
n=double(n);
nr=h5read(file,'/rand/numpair');
nerr=h5read(file,'/rand/numpair_err');
nmock=h5read(file,'/rand/numrand');
nmock=double(nmock);
rv=comoving_200b_radius(m,0.3)/1000;
% errorbar(r,n./r.^2,sqrt(n)./r.^2,[colors(i),'-']);hold on;
% plot(r,nr./r.^2,[colors(i),'--']);
rat=n./(nr.*nmock/max(nmock));
% rat_err=sqrt((1./n+(nerr./nr).^2)).*rat;
rat_err=nerr./nr;
errorbar(r/rv,rat-1,rat_err,[colors(i),'o']);
hold on;
% plot(r/rv,gammaz(z,dz).*nfw_surf_overdensity(r/(1+z),m,z,2),[colors(i),'-'])
plot(r/rv,gammaz(z,dz).*sub_surf_overdensity(r/(1+z),m,z,2,1/alphas(i)),[colors(i),'--']) %r_{-2}=0.6R200c~[0.4~0.5]R200b
% f=fittype(@(alp,x) gammaz(z,dz).*sub_surf_overdensity(x/(1+z),m,z,2,1/alp));
% fited=fit(r(r<rv&r>0.06),rat(r<rv&r>0.06)-1,f,'start',0.5,'lower',0.2,'upper',0.6)
% plot(fited,'-')
%----------------
% s1=zeros(size(r));
% s2=0;
% for gid=1:ngrps
%     if grp.Mult(gid)<=2, continue;end
%     if grp.LumMass(gid)<mmin(id)||grp.LumMass(gid)>mmax(id), continue;end
%     zh=grp.Zfof(gid);
%     mh=grp.LumMass(gid)/1e10;
% %     m=grp.MassProxy(gid)/1e10;
%     if isnan(mh)||mh==0, continue;end
% %     s1=s1+quad(@(z) normpdf(z,zh,sigz(zh)).*source_zdistr(z).*h(z)/c,zh+dz,1.5).*nfw_surf_overdensity(r/(1+zh),mh,zh,1);
%     s1=s1+quad(@(z) normpdf(z,zh,sigz(zh)).*source_zdistr(z).*h(z)/c,zh+dz,1.5);%.*sub_surf_overdensity(r/(1+zh),mh,zh,2,1/alphas(i));
%     s2=s2+quad(@(z) source_zdistr(z),zh+dz,1.5);
% end
% plot(r,s1./s2,[colors(i),':']);
% plot(r,s1./s2.*sub_surf_overdensity(r/(1+z),m,z,2,1/alphas(i)),[colors(i),':']);
% s1./s2,gammaz(z,dz)
% f=fittype(@(alp,x) s1(1)/s2.*sub_surf_overdensity(x/(1+z),m,z,2,1/alp));
% fited=fit(r(r<rv&r>0.04),rat(r<rv&r>0.04)-1,f,'start',0.5,'lower',0.3,'upper',0.6)
% plot(fited,'-')
end


set(gca,'xscale','log')
% set(gca,'yscale','log')
%%
id=2;
bid=num2str(bid);
files={['WL_L',bid,'.DZ0.hdf5'],['WL_L',bid,'.hdf5'],['WL_L',bid,'.DZ2.hdf5'],['WL_L',bid,'.DZ3.hdf5']};
dzs=[0.01,0.1,0.2,0.3];
alphas=[0.5,0.4,0.3,0.3];
y=dz;
yg=y;
for i=1:numel(files)
file=files{i}
dz=dzs(i);
% file=['WL_L4.hdf5'];
m=h5read(file,'/predict1/Mmean');
z=h5read(file,'/predict1/z');
r=h5read(file,'/shear/seperation');
n=h5read(file,'/shear/numpair');
n=double(n);
nr=h5read(file,'/rand/numpair');
nerr=h5read(file,'/rand/numpair_err');
nmock=h5read(file,'/rand/numrand');
nmock=double(nmock);
rv=comoving_200b_radius(m,0.3)/1000;
rat=n./(nr.*nmock/max(nmock));
rat_err=nerr./nr;
y(i)=sum(rat(1:1)-1);
yg(i)=gammaz(z,dz);
end
plot(dzs,y/y(1),'o');
hold on;
plot(dzs,yg/yg(1),'r-');