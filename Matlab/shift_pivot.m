function [AA,CC]=shift_pivot(par,bias,err,C,dx)
% update parameters in Table 1 of the paper
% dx=[1;log10(Pivot_new)-log10(Pivot)]. column vec
% par: [log10(normalization); slopes]. column vec
% bias: par_true-par
% err: errors on pars
% C: correlation coefficient matrix
% output: 
% AA: [A,bA,sA], new A,its bias,and error
% CC: new correlation coef between A and par

A=par'*dx;
bA=bias'*dx;
cov=err*err'.*C;
sigAP=cov*dx;
sigAA=dx'*sigAP;
sA=sqrt(sigAA);
CC=sigAP./err/sA;
CC(1)=1.;
AA=[A,bA,sA];
