c='rgbmckrgbmck';
% myfigure;
subplot(2,2,1);
zl=0.1:0.2:1.;
for i=1:numel(zl)
    zs=linspace(zl(i),2,20);
    y=1./sigma_crit(0.3,0.7,zl(i),zs);
    plot(zs,y,c(i))
    hold on;
end
xlabel('source redshift');ylabel('$1/\Sigma_{crit}$','interpreter','latex')
l=legend(mat2cell([repmat('z_l=',numel(zl),1),num2str(zl')]));
set(l,'location','northwest');
subplot(2,2,2);
zs=0.1:0.2:1.;
zmax=zs;
ymax=zs;
for i=1:numel(zs)
    zl=linspace(0.01,zs(i),50);
    y=1./sigma_crit(0.3,0.7,zl,zs(i));
    plot(zl,y,c(i))
    hold on;
    [ymax(i),j]=max(y);
    zmax(i)=zl(j);
end
xlabel('lens redshift');ylabel('$1/\Sigma_{crit}$','interpreter','latex')
l=legend(mat2cell([repmat('z_s=',numel(zs),1),num2str(zs')]));
set(l,'location','northeast');
subplot(2,2,3);
plot(zs,zmax)
hold on;
k=zmax/zs;
plot(zs,zmax/zs*zs,'r')
xlabel('source redshift')
ylabel('Most Efficient lens redshift')
l=legend('numerical',['z_m=',num2str(k,'%.1f'),'z_s']);
set(l,'location','northwest')
subplot(2,2,4);
dz=[0.1,0.3,0.5,1];
zl=0.1:0.1:1.;
for i=1:numel(dz)
    y=1./sigma_crit(0.3,0.7,zl,zl+dz(i));
    plot(zl,y,c(i))
    hold on;
end
xlabel('lens redshift');ylabel('$1/\Sigma_{crit}$','interpreter','latex')
l=legend(mat2cell([repmat('dz=',numel(dz),1),num2str(dz')]));
set(l,'location','northeast');
print('-depsc','/work/Projects/Lensing/outputv4/Theory/lensing_efficiency.eps');