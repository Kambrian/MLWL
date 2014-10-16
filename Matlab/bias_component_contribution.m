r=logspace(-2,2);
z=0.2;
color='rgbk';
i=1;
myfigure;
for M=[1e2,1e3,1e4,1e5]
    [y,ynfw,ylin,rv]=lensing_rfunc_autobias(r,M,z,-1);
    plot(r/rv,ylin./y,'-','color',color(i));
    i=i+1;
    hold on;
end
%%
set(gca,'xscale','log');
set(gca,'yscale','log');
xlim([0.1,20]);
ylim([1e-3,3]);
xlabel('$R/R_{200b}$');
ylabel('$\Delta\Sigma_{lin}/\Delta\Sigma_{tot}$');
title('z=0.2');
l=legend('$10^{12}$','$10^{13}$','$10^{14}$','$10^{15}$');
set(l,'interpreter','latex','location','northwest')
grid on;
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/biascomponent_contribution.eps');
