function [DynMass,C1,C2,A1,A2]=dynamical_mass(cat,mockmass)
% calibrate cat.MassProxy against input mass to be used as DynMass
%


% lum=[0,10:0.5:12,13];[~,ibin]=histc(log10(cat.MassProxy),lum);
mult=[2,5,10,20,1000];[~,ibin]=histc(cat.Mult,mult);
z=[0,0.1,0.2,0.3,1];[~,jbin]=histc(cat.Zfof,z);

DynMass=cat.MassProxy;
C1=zeros(4,4);
x=C1;y=C1;
C2=C1;
for i=1:4
    for j=1:4
    f=(ibin==i&jbin==j);
    f=f&mockmass>0;
%     f=f&cat.Mult>2;
    disp(sum(f))
    C1(i,j)=median(mockmass(f)./cat.MassProxy(f));
    mbin=median(mockmass(f));
    fbin=median(cat.MassProxy(f));
    C2(i,j)=mbin/fbin;
    DynMass(f)=cat.MassProxy(f)*C1(i,j);
    x(i,j)=median(1./sqrt(cat.Mult(f)));
    y(i,j)=median(1./sqrt(cat.Zfof(f)));
    end
end
X=[ones(numel(x),1),x(:),y(:)];
A1=X\C1(:);
A2=X\C2(:);
%~~Note: It's no big deal whether to fit this as a function of N and Z; rather adopt the calibration for discrete bins.
% DynMass=(A(1)+A(2)./sqrt(cat.Mult)+A(3)./sqrt(cat.Zfof)).*cat.MassProxy;

% if sum(ibin==numel(lum)|ibin==0)>0
if sum(ibin==numel(mult)|ibin==0)>0    
    error('unexpected multiplicity: >1000??');
end
if sum(jbin==numel(z)|jbin==0)>0
    error('unexpected redshift: not [0,1]?');
end
