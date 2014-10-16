function D=growth_factor(OmegaM,z)
% linear growth factor normalized to z=0, for flat cosmology
% input: OmegaM, z=0 matter density
%            z: redshift


x0=(1/OmegaM-1)^(1/3);
x=x0./(1+z);
OmegaMZ=1./(1+x.^3);
D3=z;
for i=1:numel(z)
D3(i)=(1./OmegaMZ(i)-1).^(1/3)*hypergeom([1/3,1],11/6,1-1./OmegaMZ(i));  % Nakamura & Suto 1997, eq.[C-25]
end
D0=(1./OmegaM-1).^(1/3)*hypergeom([1/3,1],11/6,1-1./OmegaM);
D=D3/D0;
