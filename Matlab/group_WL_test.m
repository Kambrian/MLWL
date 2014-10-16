function [et,ex,es,sigmc]=group_WL_test(cen,rrange,zcen,sky,flag_comoving)
% lensing output for each group
% input: cen,zcen: center [ra,dec] and redshift
%        rrange: comoving radial range to calculate the the signal
%        sky: sky patch id to search for source
%        nbin: number of radial bins
% output:
%        m,n: sum and number for M_\zeta, with M_zeta=m./n;
%        s,w: sum and sum of weight for DeltaSigma, with DeltaSigma=s./w;
%        em2,es2: sum of square errors 
%        output seperately so that you can further sum m and n (or s and w)
%        with signal from other groups before making the average
global OmegaM
switch flag_comoving  %this only affects the surface density estimator, but not the mass estimator
    case 0
        scale=1;
    case 1
        scale=1./(1+zcen);
    otherwise
        error('flag_comoving must be 0 or 1');
end

dl=AD_dist(OmegaM,0,zcen);  %angular diameter distance for lens
thetamax=rrange(2)./dl*180/pi.*scale;
thetamin=rrange(1)./dl*180/pi.*scale;
    
[et,ex,es,theta,sigmc]=get_sources(cen,thetamax,zcen,sky);%note sigmc is physical everywhere

f=theta<thetamin;
et=et(f);
ex=ex(f);
es=es(f);
sigmc=sigmc(f);