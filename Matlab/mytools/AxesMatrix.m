function ax = AxesMatrix(n,ratioh,ratiow,margin)
%  AxesMatrix: create a matrix of close axes
%                       axes are indexed from top to bottom, then left to right
%  n: MatrixSize [nh*nw]
%  ratioh: relative heights [h1,h2,h3...]
%  ratiow: relative width  [w1,w2,w3...], set to 1 to use equal size
%  margin: [l,r,b,t]  margins to the left, right, bottom, top; optional.
%
% Problem: axis location is overwritten by subsequent plot calls. why???

nh=n(1);
nw=n(2);
if nargin<4 
    margin=[0.15,0.15,0.1,0.1];
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

ax=ones(n);
for i=1:nh
    for j=1:nw
        ax(i,j)=axes('position',[xoffset(j),yoffset(i),ratiow(j),ratioh(i)]);
        if i==nh
            set(gca,'xaxislocation','bottom');
        end
        if i==1
            set(gca,'xaxislocation','top');
        end
        if j==nw
            set(gca,'yaxislocation','right');
        end
        if j==1
            set(gca,'yaxislocation','left');
        end
        if i>1&&i<nh
           set(gca,'xticklabel',[]);
        end
        if j>1&&j<nw
            set(gca,'yticklabel',[]);
        end
    end
end
set(ax,'box','on');

end

