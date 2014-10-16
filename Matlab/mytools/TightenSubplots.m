function TightenSubplots(ax,ratioh,ratiow,margin)
%  AxesMatrix: create a matrix of close axes
%                       axes are indexed from top to bottom, then left to right
%  ax: axes handels to format, has dimension [nh*nw]
%  ratioh: relative heights [h1,h2,h3...], dimension nh*1; optional
%  ratiow: relative width  [w1,w2,w3...], set to 1 to use equal size;
%               optional;
%  margin: [l,r,t,b]  margins to the left, right, bottom, top; optional.
%
% examples: TightenSubplots(ax)
%                    TightenSubplots(ax,rh)
%                    TightenSubplots(ax,rh,rw)
%                    TightenSubplots(ax,rh,rw,margin)

nh=size(ax,1);
nw=size(ax,2);
if nargin<4 
    margin=[0.15,0.15,0.1,0.1];
end

if nargin<3
    ratiow=1;
end

if nargin<2
    ratioh=1;
end

if numel(ratioh)==1
    ratioh=ones(nh,1);
end
if numel(ratiow)==1
    ratiow=ones(nw,1);
end

ratioh=ratioh/sum(ratioh)*(1-margin(3)-margin(4));  %heights
ratiow=ratiow/sum(ratiow)*(1-margin(1)-margin(2)); %widths
xoffset=cumsum([margin(1);ratiow(:)]);
xoffset=xoffset(1:end-1);
yoffset=cumsum([margin(3);ratioh(:)]);
yoffset=1-yoffset(2:end);

for i=1:nh
    for j=1:nw
        set(ax(i,j),'position',[xoffset(j),yoffset(i),ratiow(j),ratioh(i)]);
        if i==nh
            set(ax(i,j),'xaxislocation','bottom');
        end
        if i==1&&i<nh
            set(ax(i,j),'xaxislocation','top');
        end
        if j==nw
            set(ax(i,j),'yaxislocation','right');
        end
        if j==1
            set(ax(i,j),'yaxislocation','left');
        end
        if i>1&&i<nh
           set(ax(i,j),'xticklabel',[]);
           set(get(ax(i,j),'xlabel'),'string',[]);
        end
        if j>1&&j<nw
            set(ax(i,j),'yticklabel',[]);
            set(get(ax(i,j),'ylabel'),'string',[]);
        end
    end
end
set(ax,'box','on');

end

