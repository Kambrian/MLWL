function sig=sigma_crit_eff(zl)
%effective sigma_crit for S/N calculation, for lens zl
%sigc_eff=<sigc^-2>^(-1/2);

f=@(z) sigma_crit(0.3,0.7,zl,z).^-2.*source_zdistr(z);
nf=@(z) f(z)/f(zl+0.5);
% g=@(z) source_zdistr(z);
sig=quad(nf,zl+0.1,20)*f(zl+0.5);
sig=1/sqrt(sig);
