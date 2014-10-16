function s=lin_corr_DeltSig(rp)
%DeltaSigma(r), WL density profile from linear correlation function
%rp: comoving, Mpc/h

global proj_spline_pp
if isempty(proj_spline_pp)
    disp('interp proj_spline_pp...');
er=linspace(-2,4,50);
r=10.^er;
kl=proj_corr_lin(r);
% ekl=log10(kl);
% pp=spline(r,kl);
proj_spline_pp=pchip(r,kl);
end

s=rp;
for i=1:numel(rp)
    s(i)=quadl(@(x) ppval(proj_spline_pp,x).*x*2,0,rp(i))/rp(i)^2-ppval(proj_spline_pp,rp(i));
end