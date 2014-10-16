function y=signlog(x)
%log10 with sign
y=(log10(abs(x)));
y(y<0)=0;
y=sign(x).*y;