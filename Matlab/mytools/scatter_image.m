function h=scatter_image(x,y,z,type)
% to produce color-coded 2-D view from general 3-D points

x=x(:);
y=y(:);
z=z(:);
if nargin<4
    type='scatter';
end
switch type
    case 'scatter'
    ncolor=10;
    c=colormap(colorcube(ncolor));
    [~,ic]=histc(z,linspace(min(z)-eps,max(z)+eps,ncolor));
    for i=1:numel(z)
    plot(x(i),y(i),'.','color',c(ic(i),:));
    hold on;
    end
    case 'image'
    ngrid=50;
    xi=linspace(min(x),max(x),ngrid);
    yi=linspace(min(y),max(y),ngrid);
    F=TriScatteredInterp(x,y,z);
    [xx,yy]=meshgrid(xi,yi);
    zz=F(xx,yy);
    %  [xx,yy,zz]=griddata(x,y,z,xi,yi');
    imagesc(minmax(xi),minmax(yi),zz);
    set(gca,'ydir','normal');
    hold on;
    contour(xx,yy,zz);
    case 'contour'
    ngrid=50;
    xi=linspace(min(x),max(x),ngrid);
    yi=linspace(min(y),max(y),ngrid);
    F=TriScatteredInterp(x,y,z);
    [xx,yy]=meshgrid(xi,yi);
    zz=F(xx,yy);
    %  [xx,yy,zz]=griddata(x,y,z,xi,yi');
    contour(xx,yy,zz);
end