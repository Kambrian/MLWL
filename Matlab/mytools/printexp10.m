function s=printexp10(x,fmt)
% print x as a*10^b
% precision of a specified in fmt; 
% if no fmt specified, default to '%1.0f' 

if nargin<2
    fmt='%1.0f';
end

j=floor(log10(x));
i=x/10^j;
if i==1
    s=['10^{',num2str(j),'}'];
else
s=[num2str(i,fmt),'\times 10^{',num2str(j),'}'];
end