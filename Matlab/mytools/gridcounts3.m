function [xx,yy,count]=gridcounts(x,y)
% xx,yy,count : size(s(2)xs(1))
% return the counts on cells centered at (xx,yy)

xmin=min(x);
xmax=max(x);
nx=xmax-xmin+1;
ymin=min(y);
ymax=max(y);
ny=ymax-ymin+1;
ix=x-xmin+1;
iy=y-ymin+1;
[xx,yy]=meshgrid(xmin:xmax,ymin:ymax);
count=zeros(size(xx));
for i=1:numel(x)
    count(iy(i),ix(i))=count(iy(i),ix(i))+1;
end
