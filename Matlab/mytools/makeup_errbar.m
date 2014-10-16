function makeup_errbar(x,y,ey,linespec)
% to add vertical lines extending from y+ey to 0

f=y<ey;
x=x(f);
y=y(f);
ey=ey(f);
hold on;
ym=get(gca,'ylim');
for i=1:numel(x)
    plot([x(i),x(i)],[y(i)+ey(i),ym(1)],linespec);
end