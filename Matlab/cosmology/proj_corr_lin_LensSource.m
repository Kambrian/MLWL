function ksi=proj_corr_lin_LensSource(rp,zl,zs)
%projected correlation function taking into account the lensing kernel with
%lens at zl and source at zs

OmegaM=0.3;h=0.7;fb=0.168;sig8=0.8;z=0;
OmegaL=1-OmegaM;
%%
global corr_spline_pp
if isempty(corr_spline_pp)
    disp('interp corr_spline_pp...');
er=linspace(-2,5,50);
r=10.^er;
kl=matter_corr_lin(OmegaM,h,fb,sig8,r,z);
% ekl=log10(kl);
% pp=spline(r,kl);
corr_spline_pp=pchip(r,kl);
% loglog(r,kl,'.');
% hold on;
% r2=logspace(-2,5,1000);
% plot(r2,ppval(pp,r2),'-');
end
%%
% rp=logspace(-1,3,50);

tol=1e-2;
abstol=1e-6;
maxinter=1e5;

ksi=rp;
sig_crit_zl=sigma_crit(OmegaM,OmegaL,zl,zs);
dDdz=@(z) 3000./sqrt((1+z).^2.*(1+OmegaM*z)-z.*(2+z)*OmegaL);  %c/H0=3000Mpc/h
dl=comoving_dist(OmegaM,OmegaL,zl);
for i=1:numel(rp)
%     f=@(r) 2.*ppval(corr_spline_pp,r).*r./sqrt(r.^2-rp(i)^2);
%     ksi(i)=quadgk(f,rp(i),1e5,'RelTol',tol,'AbsTol',abstol,'MaxIntervalCount',maxinter);
    f=@(z) ppval(corr_spline_pp,sqrt(rp(i)^2+(comoving_dist(OmegaM,OmegaL,z)-dl).^2))./sigma_crit(OmegaM,OmegaL,z,zs)*sig_crit_zl.*dDdz(z);
    ksi(i)=quadgk(f,0,zs,'RelTol',tol,'AbsTol',abstol,'MaxIntervalCount',maxinter);
end
% figure;loglog(rp,ksi)
ksi_z=@(z,rp) ppval(corr_spline_pp,sqrt(rp^2+(comoving_dist(OmegaM,OmegaL,z)-dl).^2));
kernel_z=@(z) sig_crit_zl./sigma_crit(OmegaM,OmegaL,z,zs).*dDdz(z);%lensing kernel in redshift space
figure;
z=0:0.01:0.5;
y1=ksi_z(z,5);
y2=kernel_z(z);
y1=y1/max(y1);
y2=y2/max(y2);
plot(z,y1,'r-');
hold on;
plot(z,y2,'g--');
legend('\ksi(z)','kernel');
xlabel('z');