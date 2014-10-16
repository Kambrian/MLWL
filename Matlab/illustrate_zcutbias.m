figure;
f=@(x,mu) normpdf(x,mu,2);
x=0:0.1:8;
zcut=3;
dz=2;
z1=zcut-dz;
z2=zcut+dz;
plot(x,f(x,z1));
hold on;
plot([z1,z1],[0,f(0,0)],'--');

x=0:0.1:12;
plot(x,f(x,z2));
plot([z2,z2],[0,f(0,0)],'--');
plot([zcut,zcut],[0,f(0,0)],'-');
z3=z2+dz;
plot([z3,z3],[0,f(dz,0)],':');

xp=zcut:(12-zcut)/50:12;
for i=1:50
    plot([xp(i),xp(i)],[0,f(xp(i),z1)],':');
end

xp=z3:(12-z3)/30:12;
for i=1:30
    plot([xp(i),xp(i)],[0,f(xp(i),z2)],':');
end

plot([0,0],[0,f(0,0)]);
box off;

text(-0.2,-0.12*f(0,0),'Lens');
text(z1,-0.08*f(0,0),'z_{1}');
text(zcut,-0.12*f(0,0),'z_{cut}');
text(z2,-0.08*f(0,0),'z_{2}');
text(z3,-0.08*f(0,0),'z_{3}');

% annotation('arrow',[0,1.3],[0,0]);
ylim([-0.1*f(0,0),f(0,0)]);

set(gca,'visible','off');

print('-deps','zcut_bias.eps');
