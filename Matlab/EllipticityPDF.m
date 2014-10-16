load shearmap
e=[e{1};e{2};e{3}];
[m,s]=normfit(e(:,3));
myfigure;
x=linhist(e(:,3),100,'stairs','r-');
hold on;
dx=mean(diff(x));
plot(x,normpdf(x,m,s)*0.0403*size(e,1),'g-');
% [m,s]=normfit(e(:,4));
% x=linhist(e(:,4),100,'stairs','g-');
% dx=mean(diff(x));
% plot(x,normpdf(x,m,s)*0.0403*size(e,1),'g-');

xlabel('$e_1$');ylabel('Counts');
print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/extra/ellipticityPDF.eps');
%%
[m,s]=normfit(e(:,3));
myfigure;
histfit(sqrt(e(:,3).^2+e(:,4).^2),100,'rayleigh')
% x=linhist(e(:,3),100,'stairs','r-');
hold on;
% dx=mean(diff(x));
% plot(x,normpdf(x,m,s)*0.0403*size(e,1),'g-');
% [m,s]=normfit(e(:,4));
% x=linhist(e(:,4),100,'stairs','g-');
% dx=mean(diff(x));
% plot(x,normpdf(x,m,s)*0.0403*size(e,1),'g-');

xlabel('$e$');ylabel('Counts');
% print('-depsc','/work/Projects/Lensing/outputv4/paper/v3up8/extra/ellipticityPDF.eps');