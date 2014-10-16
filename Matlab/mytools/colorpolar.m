function colorpolar(t,r,c)

t=t/180*pi;

% 
polar(0,max(r)*1.1);
set(gco,'visible','off');
hold on;
for i=1:numel(t)
    plot(r(i).*cos(t(i)),r(i).*sin(t(i)),'o','markersize',8,'markeredgecolor',c(i,:));
    hold on;
end
hold off;