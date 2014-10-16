function [w,ew,nmax]=weight_average_correction(w,ew,n)
% correct for averaging-bias
nmax=max(n);
ew=(ew.^2+w.^2).*n/nmax;  % <w^2>
w=w.*n/nmax;
ew=sqrt(ew-w.^2);  