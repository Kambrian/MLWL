function [TS,p]=Sigma2TS(sig,dof)
% TS,p=Sigma2TS(sig,dof)
% return the TS corresponding to sig-sigma gaussian significance
% assume DoF=1 by default

if nargin<2
    dof=1;
end

p=2*normcdf(sig,0,1)-1;

TS=chi2inv(p,dof);