function draw_rect(xminmax,yminmax,linespec)

x=xminmax;y=yminmax;
plot([x(1),x(1),x(2),x(2),x(1)],[y(1),y(2),y(2),y(1),y(1)],linespec);