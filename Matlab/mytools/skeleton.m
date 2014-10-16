function [xmed,ymed,ylim,xm,ym,ysig,count]=skeleton(x,y,nbin,alpha)
% [xmed,ymed,ylim,xm,ym,ysig,count]=skeleton(x,y,nbin,alpha)
% to divide x into bins and give estimation of center and variance of y
% inside each bin
%input:
% x,y: column vectors to extract skeleton from
% nbin: number of bins or bin edges for x
% alpha: confidence level for boundary estimation
%

if nargin<4
    alpha=0.683;
end

if numel(nbin)>1
    xbin=nbin;
    nbin=numel(nbin);
else
xbin=linspace(min(x)-2*eps,max(x)+2*eps,nbin)';
end

[count,bin]=histc(x,xbin);

xm=zeros(nbin-1,1);
ym=xm;
ysig=xm;
xmed=xm;
ymed=xm;
ylim=zeros(nbin-1,2);
% alpha=(1-alpha)/2;
for i=1:nbin-1
    xm(i)=mean(x(bin==i));
    xmed(i)=median(x(bin==i));
    ym(i)=mean(y(bin==i));
    ymed(i)=median(y(bin==i));
    ysig(i)=std(y(bin==i));
%     tmp=sort(y(bin==i));
    if count(i)
%     ylim(i,:)=[tmp(ceil(alpha*count(i))),tmp(ceil((1-alpha)*count(i)))];
    ylim(i,:)=prctile(y(bin==i),[0.5-alpha/2,0.5+alpha/2]*100);
    else
     ylim(i,:)=[NaN,NaN];
    end
end
count(end)=[];