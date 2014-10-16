function sig=noise_eff(zl,dz, flag_comov)
%effective noise for S/N calculation of DeltaSig, for lens zl, buffer zone width dz
%Noise^2=1/sum(w_i)=1/sum(1/sigma_i^2), so
%(S/N)^2~sum(1/sigma^2)~sum(1/sigma_crit^2)~\int sigma_crit^{-2} dN/dz dz
% if stack in comoving DSig, then (S/N)^~ \int sigma_crit^{-2} scaleF^{-4}
% dN/dz dz.
% return 1/N
% flag_comov: 0, physical (default); 1, comoving Dsig.

if nargin<3
    flag_comov=0;
end
% f=@(z) sigma_crit(0.3,0.7,zl,z).^-2.*source_zdistr(z);
% g=@(z) source_zdistr(z);
zmin=zl+dz;
if numel(zmin)>1
if numel(zl)==1
    zl=repmat(zl,size(dz));
end
if numel(dz)==1
    dz=repmat(dz,size(zl));
end
end
zmax=1.5;
sig=zeros(size(zmin));
for i=1:numel(zmin)
    if zmin(i)<zmax
        if flag_comov
        sig(i)=quad(@(z) sigma_crit(0.3,0.7,zl(i),z).^-2.*source_zdistr(z).*(1+zl(i)).^4,zmin(i),zmax);
        else
        sig(i)=quad(@(z) sigma_crit(0.3,0.7,zl(i),z).^-2.*source_zdistr(z),zmin(i),zmax);
        end
    end
end
sig=sqrt(sig);
