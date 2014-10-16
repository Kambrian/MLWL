x=1:5:50;
dof=3;
colors=colormap(jet(numel(x)));
p=TS2Sigma(x,dof);
% plot(x,p,'--');
% hold on;
dx=0.1:0.2:5;
for i=1:numel(x)
plot(dx,TS2Sigma(x(i)+dx,dof+1)/p(i),'color',colors(i,:)); hold on;
end
%%
dof=1:6;
colors=colormap(jet(numel(dof)));
for j=1:numel(dof)
    x=dof(j):5:50;
    dx=x;
    for i=1:numel(x)
        dx(i)=fzero(@(dx) chi2cdf(x(i)+dx,dof(j)+1)/chi2cdf(x(i),dof(j))-1,1);
    end
    plot(x,dx,'color',colors(j,:));hold on;
end
xlabel('TS');
ylabel('dTS');
