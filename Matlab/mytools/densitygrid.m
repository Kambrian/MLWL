function [xx,yy,count,cellsize]=densitygrid(x,y,s,xrange,yrange,w)
% xx,yy,count : size(s(2)xs(1))
% return the counts on cells centered at (xx,yy)
% xrange, yrange optional.
%
if nargin==5
    xmin=xrange(1);
    xmax=xrange(2);
    ymin=yrange(1);
    ymax=yrange(2);
else
xmin=min(x);
xmax=max(x);
ymin=min(y);
ymax=max(y);
end
if nargin<6
    w=ones(size(x));
end
f=x>xmin&x<xmax&y>ymin&y<ymax;x=x(f);y=y(f); %trim data
hx=(xmax-xmin)/s(1);
hy=(ymax-ymin)/s(2);
ix=ceil((x-xmin)/hx);
iy=ceil((y-ymin)/hy);
% ix(ix<=0)=1;ix(ix>s(1))=s(1);
% iy(iy<=0)=1;iy(iy>s(2))=s(2);
[xx,yy]=meshgrid(xmin+(0.5:s(1)-0.5)*hx,ymin+(0.5:s(2)-0.5)*hy);
iz=sub2ind([s(2),s(1)],iy,ix);
[count,bin]=histc(iz,(0:s(1)*s(2))+0.5);
for i=1:numel(count)
    count(i)=sum(w(bin==i));
end
count(end)=[];
count=reshape(count,[s(2),s(1)]);
cellsize=hx*hy;
