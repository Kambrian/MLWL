sigz=@(z) 1.6./(z.^2-14.5*z+21.7);
H=100;
OmegaM=0.3;
h=@(z) H*sqrt(OmegaM .*(1+z).^3+(1-OmegaM));
c=3e5;
G=43007.1;
gammaz=@(zh,dz) quad(@(z) normpdf(z,zh,sigz(zh)).*source_zdistr(z).*h(z)/c,zh+dz,1.5)./quad(@(z) source_zdistr(z),zh+dz,1.5);
%%

dz=0:0.05:0.5;
zl=0.2;
zref=0.01;
sn=noise_eff(zl,dz)/noise_eff(zl,zref);
plot(dz,sn,'r')
hold on;
y=dz;
for i=1:numel(dz)
    y(i)=gammaz(zl,dz(i));
end
plot(dz,y/gammaz(zl,zref),'g');
legend('SN','Contamination');
% set(gca,'yscale','log')
xlabel('buffer redshift width');
ylabel('SN or Contamination Amplitude');
print('-depsc','/work/Projects/Lensing/outputv4/zbuffer_optimize.eps')
%%
dz=0:0.05:0.5;
zl=0:0.05:0.5;
zref=0.3;
c='rgbk'
plot(zl,noise_eff(zl,zref,0),'r')
hold on;
plot(zl,noise_eff(zl,zref,1),'g')
%%
for i=1:numel(zl)
sn=noise_eff(zl(i),dz);
plot(dz,sn,c(i))
hold on;
end
