function [HybMass,C]=luminosity_mass(cat)
% calibrate TotFluxProxy to be used as LumMass
%


mult=[2,5,10,20,1000];
z=[0,0.1,0.2,0.3,0.5];
[~,ibin]=histc(cat.Nfof,mult);
[~,jbin]=histc(cat.Zfof,z);

HybMass=cat.HybMass;
C=zeros(4,4);
for i=1:4
    for j=1:4
    f=(ibin==i&jbin==j);
    disp(sum(f))
    mbin=median(cat.MassProxy(f));
    fbin=median(cat.HybMass(f));
    C(i)=mbin/fbin;
    HybMass(f)=cat.HybMass(f)*C(i);
    end
end

if sum(ibin==5|ibin==0)>0
    error('unexpected multiplicity: >1000??');
end
if sum(jbin==5|jbin==0)>0
    error('unexpected redshift: not [0,0.5]?');
end
