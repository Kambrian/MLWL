function [xm,ym,dyx,h,xr]=linhist(data,nbin,plottype,linespec)
% [xm,ym,dyx,h,xr]=linhist(data,nbin,plottype,linespec)
% plottype: 'stairs','stairsnorm','norm', or empty ''for plain plot.
if numel(nbin)>1
    x=nbin;
    nbin=numel(nbin);
else
x=linspace(min(data)-2*eps,max(data)+2*eps,nbin)';
end
xr=x(2:end);

[count,bin]=histc(data,x);

xm=x(1:nbin-1);
for i=1:nbin-1
    xm(i)=mean(data(bin==i));
end
% xm(isnan(xm))=x([isnan(xm);false]);

count=reshape(count,size(x));
ym=count(1:nbin-1);
dyx=ym./diff(x);
h=0;
if nargin>2
    if nargin==3,linespec='-'; end
    if strcmp(plottype,'stairs')
        h=stairs(x,count,linespec);
    elseif strcmp(plottype,'stairsnorm')
        h=stairs(x,count/sum(count),linespec);
    elseif strcmp(plottype,'norm')
        h=plot(xm,ym/sum(count),linespec);
    elseif isempty(plottype)||strcmp(plottype,' ')
        h=plot(xm,ym,linespec);
    else
        eval(['h=',plottype,'(xm,ym)']);
    end
end
end

