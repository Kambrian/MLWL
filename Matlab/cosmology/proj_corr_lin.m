function ksi=proj_corr_lin(rp)

OmegaM=0.3;h=0.7;fb=0.168;sig8=0.8;z=0;
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
for i=1:numel(rp)
    f=@(r) 2.*ppval(corr_spline_pp,r).*r./sqrt(r.^2-rp(i)^2);
    ksi(i)=quadgk(f,rp(i),1e5,'RelTol',tol,'AbsTol',abstol,'MaxIntervalCount',maxinter);
end
% figure;loglog(rp,ksi)