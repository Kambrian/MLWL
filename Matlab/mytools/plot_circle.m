function h=plot_circle(cen,r)

t=0:0.1:2*pi+0.1;
x=cen(1)+r*sin(t);
y=cen(2)+r*cos(t);
h=line(x,y);