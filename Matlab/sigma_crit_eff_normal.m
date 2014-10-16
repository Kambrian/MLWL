function sig=sigma_crit_eff_normal(zl)
%effective sigma_crit for S/N calculation, for lens zl
% sigc_eff= sqrt(<sigc^2>)/fraction(z)

f=@(z) sigma_crit(0.3,0.7,zl,z).^2.*source_zdistr(z);
nf=@(z) f(z)/f(zl+0.5);
g=@(z) source_zdistr(z);
sig=quad(nf,zl+0.1,20)*f(zl+0.5);
sig=sig/quad(g,zl+0.1,20)^2;
sig=sqrt(sig);
