function p=source_zdistr(z)
%redshift distribution for source galaxies from Madalbaum 2005
% 
alpha=4;
zs=0.1; 
% alpha=2; zs=0.12 %seems to fit sdss's z
% alpha=4;zs=0.05 %seems to fit gama's z


p=z.^(alpha-1).*exp(-z./zs)./zs^alpha./gamma(alpha);

%redshift distribution for source galaxies from Madalbaum 2008

% alpha=2.35;
% zs=0.3;
% p0=(sqrt(2)*zs)^alpha/2*gamma(alpha/2);
% p=z.^(alpha-1).*exp(-0.5*(z/zs).^2)/p0;
