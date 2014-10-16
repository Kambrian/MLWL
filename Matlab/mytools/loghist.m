function [xm,ym,dyx,x]=loghist(data,nbin,plottype,linespec,weight, dim)
% function [xm,ym,dyx]=loghist(data,nbin,plottype,linespec,weight, dim)
%  weight: weight for y counts
%  dim: dyx=ym./diff(x.^dim), density in dim-dimension space
if numel(nbin)==1
x=logspace(log10(min(data(data>0))),log10(max(data)*1.01),nbin)';
% x=linspace((min(data(logical(data)))),(max(data)*1.01),nbin)';
else
x=nbin;
nbin=numel(x);
end

if nargin<5||numel(weight)<2
    weight=ones(size(data));
end

if nargin<6
    dim=1;
end

[count,bin]=histc(data,x);

xm=x(1:nbin-1);
for i=1:nbin-1
    xm(i)=mean(data(bin==i));
    count(i)=sum(weight(bin==i));
end
count(end)=sum(weight(bin==nbin));

count=reshape(count,size(x));
ym=count(1:nbin-1);
% dyx=ym./diff(x.^dim);
dyx=ym./diff(log(x));
if nargin>2
    if nargin==3||isempty(linespec),linespec='-'; end
    if strcmp(plottype,'stairs')
        stairs(x,count,linespec);
        set(gca,'xscale','log');
        set(gca,'yscale','log');
    elseif strcmp(plottype,'norm')
        loglog(xm,ym/sum(count),linespec);
    elseif isempty(plottype)||strcmp(plottype,' ')
              %do nothing
    else
    loglog(xm,ym,linespec);
    end
end

