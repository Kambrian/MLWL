data=importdata('/work/Projects/Lensing/code/v8.0/Matlab/cosmology/dndlogM_ST_3778_extra.dat',' ',1);
hmf=data.data;
%%
G=43.0071;
H0=100;
rhoc=3*H0^2/8/pi/G;
rhoc=rhoc*1e10; %Msun/h unit
omegaC=0.3*0.85;
rhom=rhoc*omegaC;
%%
figure;
plot(hmf(:,1),hmf(:,2),'o-');
% yscale('log');
xscale('log');
figure;
plot(hmf(:,1),hmf(:,7)/rhom,'o-');
% yscale('log');
xscale('log');
figure;
plot(hmf(:,1),hmf(:,6)/rhom,'o-');
% yscale('log');
xscale('log');
%%
spm=spline(log10(hmf(:,1)),log(hmf(:,2)));
sp=spline(log10(hmf(:,1)),log(hmf(:,3)));
halomfunc=@(logm) exp(ppval(sp,logm)); %dn/dlog10(M)/dV, for halo mass func
mhalomfunc=@(logm) exp(ppval(spm,logm));%m*dn/dlog10(m)/dV
%%
spm=pchip(log10(hmf(:,1)),hmf(:,2));
spcum=spline(log10(hmf(:,1)),hmf(:,7));
mhalomfunc=@(logm) ppval(spm,logm);%m*dn/dlog10(m)/dV
mcumfunc=@(logm) ppval(spcum,logm); %m(>m)
%%
figure;
loglog(hmf(:,1),hmf(:,2),'o');
hold on;
m=logspace(0,18);
loglog(m,mhalomfunc(log10(m)),'kx-')
%%
m=logspace(9,16,10);
for i=1:numel(m)
    y(i)=quad(mhalomfunc,3,log10(m(i)));
end
figure;
loglog(hmf(:,1),hmf(:,7)/rhoc);
hold on;
loglog(m,y/rhoc,'o');
%%
rhobar=quadgk(mhalomfunc,3,16)/1e10

rhobar/rhoc
%%
z=0;
m=logspace(-5,16);
b=linear_bias(m/1e10,z,0.8);
spb=spline(log10(m),b);
bias_func=@(logm) ppval(spb, logm);
%%
figure;
loglog(m,b)
yscale('linear');
%% the extra bias factor if only halos <M are contributing to the correlation function.
m=logspace(12,16);
y=m;
for i=1:numel(m)
    y(i)=quadgk(@(logm) mhalomfunc(logm).*bias_func(logm),3,log10(m(i)));
end
myfigure;
loglog(m,y/y(end),'r-');
hold on;
% loglog(hmf(:,1),hmf(:,7)/rhom,'g-');
yscale('linear');
xlabel('$M[M_\odot/h]$');
ylabel('$b_2(<M)$');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/extra/missedbias.eps');