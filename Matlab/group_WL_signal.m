function [rs,er2,m,n,s,w,em2,es2, s_pred]=group_WL_signal(cen,rrange,zcen,sky,nbin,flag_comoving,mproxy)
% lensing output for each group
% input: cen,zcen: center [ra,dec] and redshift
%        rrange: radial range to calculate the the signal, comoving or
%                                         physical depending on the flag,
%                                         Mpc/h unit.
%        sky: sky patch id to search for source
%        nbin: number of radial bins
% output:
%        rs: weighted sum of radial seperations, comoving or physical depending on
%        flag
%        er2: wighted sum of r^2, for error anal
%        m,n: sum and number for M_\zeta, with M_zeta=m./n;
%        s,w: weighted sum and sum of weight for DeltaSigma, with DeltaSigma=s./w;
%        em2,es2: sum of square errors 
%        output seperately so that you can further sum m and n (or s and w)
%        with signal from other groups before making the average
%
% mproxy: mass proxy
% s_pred: predicted s according to mproxy.
global OmegaM
switch flag_comoving  %this only affects the surface density estimator, but not the mass estimator
    case 0
        scale=1;
    case 1
        scale=1./(1+zcen);
    otherwise
        error('flag_comoving must be 0 or 1');
end

dl=AD_dist_flat(OmegaM,0,zcen);  %angular diameter distance for lens
thetamax=rrange(2)./dl*180/pi.*scale;
thetamin=rrange(1)./dl*180/pi.*scale;
if thetamax>90, thetamax=90;end % we do not need to search half of the sky

[et,ex,es,theta,sigmc,ed2]=get_sources(cen,thetamax,zcen,sky);%note sigmc is physical everywhere
% wght=es.^-2.*sigmc.^-2;
wght=es.^-2.*sigmc.^-2*scale^-4; %Mandelbaum's weighting
% wght=ones(size(et));
% wght=es.^-2;
% wght=sigmc.^-2;
% tbin=linspace(thetamin,thetamax,nbin+1);
tbin=logspace(log10(thetamin),log10(thetamax),nbin+1);
[n,bin]=histc(theta,tbin);
n(end)=[];

m=zeros(nbin,1);s=m;w=m;em2=m;es2=m;rs=m;er2=m;
s_pred=s;

if numel(theta)<=1, n=reshape(n,nbin,1);end
theta=theta/thetamax;

dsig=nfw_DeltSig(theta*rrange(2)*scale,mproxy,zcen,2,0);

for i=1:nbin
    if n(i)
        rs(i)=sum(theta(bin==i).*wght(bin==i));
        er2(i)=sum(theta(bin==i).^2.*wght(bin==i));
        tmp=et(bin==i).*sigmc(bin==i)./theta(bin==i).^2; 
        m(i)=sum(tmp);
        em2(i)=sum(tmp.^2.*ed2(bin==i));
        tmp=et(bin==i).*wght(bin==i).*sigmc(bin==i);
        s(i)=sum(tmp);
        es2(i)=sum(tmp.^2.*ed2(bin==i));
        w(i)=sum(wght(bin==i));
        tmp=dsig(bin==i).*wght(bin==i);
        s_pred(i)=sum(tmp);
    end
end
rs=rs*rrange(2);
er2=er2*rrange(2)^2;
s=s*scale^2;  %translating from physical to comoving if required by flag
s_pred=s_pred*scale^2;
% es2=w;
es2=es2*scale^4;
m(end:-1:1)=cumsum(m(end:-1:1));
em2(end:-1:1)=cumsum(em2(end:-1:1));
km=pi.*(tbin(1:end-1)'*pi/180*dl).^2; %coefficients for each halo
m=m.*km;
em2=em2.*km.^2;
n(end:-1:1)=cumsum(n(end:-1:1));

%======== To add:
% shear responsivity
% other estimators
