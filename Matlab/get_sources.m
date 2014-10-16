function [et,ex,es,theta,sigmc,ed2]=get_sources(cen,thetamax,zcen,sky)

global ll e OmegaM OmegaL  macros

sources=select_grids(cen,thetamax,[ll{sky}.xrange(1),ll{sky}.yrange(1)],ll{sky}.step,ll{sky}.ngrid,ll{sky}.grids);
% sources(e{sky}(sources,6)<zcen+0.1|e{sky}(sources,7)<zcen)=[];  % exclude close pairs
sources(e{sky}(sources,6)<zcen+0.1)=[];  
[et,ex,theta,f]=txshear(cen,[e{sky}(sources,1),e{sky}(sources,2)],[e{sky}(sources,3),e{sky}(sources,4)],thetamax);
sources(f)=[];
sigmc=sigma_crit(OmegaM,OmegaL,zcen,e{sky}(sources,6));
es=e{sky}(sources,5).^2+0.4.^2;%add intrinsic noise;  
ed2=es./et.^2;

if macros.PHOTOZ_DENS_ERR2_MAX>0
sigmcerr2=-[sigmc-sigma_crit(OmegaM,OmegaL,zcen,e{sky}(sources,7)),sigma_crit(OmegaM,OmegaL,zcen,e{sky}(sources,8))-sigmc]; %what if source-z is smaller than lens?
sigmcerr2=mean(sigmcerr2.^2,2);  %square error
eg2=sigmcerr2./sigmc.^2;

if macros.PHOTOZ_DENS_ERR_INC
ed2=ed2+eg2;  % square relative error from et and sigmc
end

if macros.PHOTOZ_DENS_ERR2_MAX~=inf
f=eg2<macros.PHOTOZ_DENS_ERR2_MAX;   %only select over redshift, never select over ellipticity, otherwise you get bias!
et=et(f);
ex=ex(f);
es=es(f);
theta=theta(f);
sigmc=sigmc(f);
ed2=ed2(f);
end
end