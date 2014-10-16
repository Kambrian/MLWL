function [LumMass,C,A]=luminosity_mass(cat)
% calibrate Luminosity to be used as LumMass
%


% lum=[0,10:0.5:12,13];[~,ibin]=histc(log10(cat.Luminosity),lum);
mult=[2,5,10,20,1000];[~,ibin]=histc(cat.Mult,mult);
z=[0,0.1,0.2,0.3,0.6];[~,jbin]=histc(cat.Zfof,z);

LumMass=cat.Luminosity;
C=zeros(4,4);
x=C;y=C;
for i=1:4
    for j=1:4
    f=(ibin==i&jbin==j);
 %   f=f&cat.Mult>2; %this filter does not necessarily make your result better, because your dynamical mass is calibrated with [2,5]; so if only calibrate M_L/M_d in [3,5], you introduce extra bias
 % but note that the (dyn-)mass is already biased(?) with Aeron's calibration and
 % selecting only Mult>2.
   % disp(sum(f))
    mbin=median(cat.DynMass(f));
    fbin=median(cat.Luminosity(f));
    C(i,j)=mbin/fbin;
%     C(i,j)=median(cat.DynMass(f)./cat.Luminosity(f));
    LumMass(f)=cat.Luminosity(f)*C(i,j);
    x(i,j)=median(1./sqrt(cat.Mult(f)));
    y(i,j)=median(1./sqrt(cat.Zfof(f)));
    end
end
X=[ones(numel(x),1),x(:),y(:)];
A=X\C(:);
%~~Note: It's no big deal whether to fit this as a function of N and Z; rather adopt the calibration for discrete bins.
% LumMass=(A(1)+A(2)./sqrt(cat.Mult)+A(3)./sqrt(cat.Zfof)).*cat.Luminosity;

% if sum(ibin==numel(lum)|ibin==0)>0
if sum(ibin==numel(mult)|ibin==0)>0    
    error('unexpected multiplicity: >1000??');
end
if sum(jbin==numel(z)|jbin==0)>0
    warning('unexpected redshift: not [0,0.5]?');
end
