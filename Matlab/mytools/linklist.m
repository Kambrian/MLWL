function [id,xrange,yrange,step]=linklist(x,y,s)

ind=(1:numel(x))';
xmin=min(x);
xmax=max(x);
hx=(xmax-xmin)/s(1);
ymin=min(y);
ymax=max(y);
hy=(ymax-ymin)/s(2);
ix=ceil((x-xmin)/hx);
iy=ceil((y-ymin)/hy);
ix(ix<=0)=1;ix(ix>s(1))=s(1);
iy(iy<=0)=1;iy(iy>s(2))=s(2);
id=cell(s);
for i=1:s(1)
    for j=1:s(2)
        id{i,j}=ind(ix==i&iy==j);
    end
end
xrange=[xmin,xmax];
yrange=[ymin,ymax];
step=[hx,hy];