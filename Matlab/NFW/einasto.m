function rho=einasto(x, alpha, xmax)
% Einasto Profile: rho=rhos*exp{-2/alpha*[(r/rs)^alpha-1]},
% with typical alpha~0.16
% input: x=r/rs
%            xmax, truncation radius, beyond which rho=0.
% output: rho/rhos
% rhos, rs usually quoted as rho_{-2}, r_{-2}.

if nargin<2
    alpha=0.16;
end

rho=exp(-2./alpha.*(x.^alpha-1));

if nargin==3
    rho(x>xmax)=0;
end