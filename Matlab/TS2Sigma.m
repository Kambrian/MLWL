function [sig,p]=TS2Sigma(TS,dof)
% p=TS2Sigma(TS,dof)
% return the probability that TS is not a fluctuation
% assume DoF=1 by default

if nargin<2
    dof=1;
end

TS(TS<0)=0;
p=chi2cdf(TS,dof);

sig=norminv(0.5+p/2,0,1);