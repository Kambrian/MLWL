function s=einasto_SurfDen(rp,alpha,rmax)
% int(rho,dz)
%where rho=einasto(r)
%           r^2=rp^2+z^2, in units of r_{-2}
% the resulting s would have the unit of rho_{-2}*rs_{-2}
% rmax: truncation radius for einasto profile, beyond which rho=0. in units
% of rs.

if nargin<2
    alpha=0.16;
end
if nargin<3
    rmax=50;
end

s=rp;
for i=1:numel(rp)
    s(i)=quadgk(@(z) einasto(sqrt(z.^2+rp(i)^2),alpha,rmax),0,rmax);
%     s(i)=quadgk(@(r) einasto(r,alpha)./sqrt(r.^2-rp(i)^2),rp(i),rmax);
end
s=s*2;