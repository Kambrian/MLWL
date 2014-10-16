function sig=sigma_crit_eff_cfht(zl)
%effective sigma_crit for S/N calculation, for lens zl
%sigc_eff=<sigc^-2>^(-1/2);

global cfht_z
sig=sum(sigma_crit(0.3,0.7,zl,cfht_z(cfht_z>zl+0.1)).^-2)/sum(cfht_z>zl+0.1);
frac=sum(cfht_z>zl+0.1)/numel(cfht_z);
sig=1/sqrt(sig*frac);
