function [xx,yy,zz,count,cellsize]=densitygrid3(x,y,z,s,xrange,yrange,zrange)
% xx,yy,count : size(s(2)xs(1)xs(3))
% return the counts on cells centered at (xx,yy)
% xrange, yrange optional.
%
if nargin==7
    xmin=xrange(1);
    xmax=xrange(2);
    ymin=yrange(1);
    ymax=yrange(2);
    zmin=zrange(1);
    zmax=zrange(2);
else
xmin=min(x);
xmax=max(x);
ymin=min(y);
ymax=max(y);
zmin=min(x);
zmax=max(y);
end
hx=(xmax-xmin)/s(1);
hy=(ymax-ymin)/s(2);
hz=(zmax-zmin)/s(3);
ix=ceil((x-xmin)/hx);
iy=ceil((y-ymin)/hy);
iz=ceil((z-zmin)/hz);
ix(ix<=0)=1;ix(ix>s(1))=s(1);
iy(iy<=0)=1;iy(iy>s(2))=s(2);
iz(iz<=0)=1;iz(iz>s(3))=s(3);
[xx,yy,zz]=meshgrid(xmin+(0.5:s(1)-0.5)*hx,ymin+(0.5:s(2)-0.5)*hy,zmin+(0.5:s(3)-0.5)*hz);
in=sub2ind([s(2),s(1),s(3)],iy,ix,iz);
count=histc(in,(0:s(1)*s(2)*s(3))+0.5);
count(end)=[];
count=reshape(count,[s(2),s(1),s(3)]);
cellsize=[hy, hx, hz];
